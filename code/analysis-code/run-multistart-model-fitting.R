# Script: run-multistart-model-fitting.R
# Purpose: explore multiple combinations of initial conditions and optimizers
#          for the QSP model fit in parallel, retain the best results, and
#          automatically generate a best-fit diagnostic figure.
#
# Usage:   Adjust the configuration block below (number of initial points,
#          optimizers to try, etc.) and run via `Rscript run-multistart-model-fitting.R`.
#
# Notes:   The script reuses the same data/model setup as run-model-fitting.R
#          but augments it with a multi-start, multi-solver workflow. All
#          results are written to `results/output/` and figures to
#          `results/figures/`.

rm(list = ls(all.names = TRUE))

# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------
library(here)
library(dplyr)
library(tidyr)
library(nloptr)
library(deSolve)
library(lhs)
library(future)
library(future.apply)
library(purrr)
library(tibble)
library(readr)

# plotting helper for the final figure
source(here::here("code/plotting-code/timeseries-plot-function.R"))

# ODE model and objective function
source(here::here("code/analysis-code/model-simulator-function.R"))
source(here::here("code/analysis-code/fit-model-function.R"))
# Dose prediction simulator (used for figure generation)
source(here::here("code/analysis-code/dose-predictions-simulator-function.R"))

# -----------------------------------------------------------------------------
# Configuration – tweak as needed
# -----------------------------------------------------------------------------
set.seed(20240529)

n_random_inits <- 20       # number of random Latin-hypercube initial guesses
use_bestfit_as_start <- TRUE
solver_candidates <- c(
  "NLOPT_LN_COBYLA",
  "NLOPT_LN_SBPLX",
  "NLOPT_LN_NELDERMEAD"
)

max_steps <- 1500
max_time_seconds <- 6 * 60 * 60
ftol_rel <- 1e-8
print_level_default <- 0

# Parallel workers (auto-clamped to available cores)
max_workers <- 24

# Settings for the ODE solver when generating predictions for figures
solvertype <- "vode"
tols <- 1e-9
dt <- 0.02
tfinal <- 7

# Doses to visualise in the final figure (match experimental scenarios)
figure_doses <- c(0, 10, 100)
figure_dose_labels <- c("no drug", "10 mg/kg", "100 mg/kg")

# -----------------------------------------------------------------------------
# Data + parameter setup (mirrors run-model-fitting.R)
# -----------------------------------------------------------------------------
fitdata <- readr::read_csv(
  here::here("data/processed-data/processeddata.csv"),
  show_col_types = FALSE
)

fixedparsdata <- readr::read_csv(
  here::here("data/processed-data/fixed-parameters.csv"),
  show_col_types = FALSE
)

fitdata <- fitdata %>%
  mutate(
    Scenario = factor(
      Scenario,
      levels = c("NoTreatment", "PanCytoVir10mg", "PanCytoVir100mg")
    ),
    Quantity = factor(
      Quantity,
      levels = c("LogVirusLoad", "IL6", "WeightLossPerc")
    ),
    Dose = c(0, 10, 100)[as.numeric(Scenario)],
    xvals = Day
  )

scenarios <- levels(fitdata$Scenario)

# Initial conditions (same as legacy script)
Y0 <- c(
  Ad = 0,
  Ac = 0,
  At = 0,
  U = 1e7,
  I = 0,
  V = 1,
  F = 0,
  A = 0,
  S = 0
)

# Baseline parameter guesses and bounds -------------------------------------------------
par_defs <- tribble(
  ~name,              ~start, ~lower,  ~upper,
  "b",                1e-8,    1e-12,   1e-2,
  "k",                1e-4,    1e-8,    1e5,
  "p",                2e3,     1e-3,    1e10,
  "kF",               1,       1e-10,   1e3,
  "cV",               10,      0.01,    1e5,
  "gF",               0.1,     1e-3,    1e3,
  "Fmax",             5,       0.1,     10,
  "hF",               1,       1e-3,    1e5,
  "gS",               10,      1e-4,    1e4,
  "cS",               1,       1e-2,    1e4,
  "Emax_F",           0.5,     1e-3,    1,
  "C50_F",            1,       1e-7,    1e6,
  "C50_V",            1,       1e-7,    1e6
)

var_by_qty <- fitdata %>%
  group_by(Quantity) %>%
  summarise(v = var(Value, na.rm = TRUE), .groups = "drop")

var_lookup <- setNames(var_by_qty$v, var_by_qty$Quantity)

sigma_all <- c(
  sigma_add_LogVirusLoad = sqrt(var_lookup[["LogVirusLoad"]]),
  sigma_prop_LogVirusLoad = 0,
  sigma_add_IL6 = sqrt(var_lookup[["IL6"]]),
  sigma_prop_IL6 = 0,
  sigma_add_WeightLossPerc = sqrt(var_lookup[["WeightLossPerc"]]),
  sigma_prop_WeightLossPerc = 0
)

sigma_to_fit <- c(
  "sigma_add_LogVirusLoad",
  "sigma_add_IL6",
  "sigma_add_WeightLossPerc"
)

sigma_fit_ini <- sigma_all[sigma_to_fit]
sigma_fixed <- sigma_all[setdiff(names(sigma_all), sigma_to_fit)]

par_ini <- c(par_defs$start, sigma_fit_ini)
names(par_ini) <- c(par_defs$name, names(sigma_fit_ini))

lower_bounds <- c(par_defs$lower, rep(1e-6, length(sigma_fit_ini)))
upper_bounds <- c(par_defs$upper, rep(1e3, length(sigma_fit_ini)))

if (any(lower_bounds <= 0)) {
  stop("All lower bounds must be positive to enable log-scale Latin hypercube sampling.")
}

# Optional: reuse previous best fit as one of the starting points
if (use_bestfit_as_start && file.exists(here::here("results", "output", "bestfit.Rds"))) {
  old_best <- readRDS(here::here("results", "output", "bestfit.Rds"))
  if (length(old_best) >= 1 && is.list(old_best[[1]]) && !is.null(old_best[[1]]$solution)) {
    par_ini <- old_best[[1]]$solution
    names(par_ini) <- old_best[[1]]$fitparnames
  }
}

fitparnames <- names(par_ini)

# Parameter labels (kept for downstream reporting)
parlabels <- c(
  "Virus infection rate",
  "Adaptive response clearance rate",
  "Virus production rate",
  "Innate response suppression strength",
  "Virus removal rate",
  "Maximum innate response induction",
  "Maximum innate response",
  "Adaptive response half-maximum induction",
  "Symptom induction rate",
  "Symptom decay rate",
  "Maximum innate response suppression",
  "Half maximum of innate response effect",
  "Half maximum of virus suppression effect",
  "Sigma of LogVirusLoad",
  "Sigma of IL6",
  "Sigma of WeightLossPerc"
)

# Fixed parameter samples (baseline + optional LHS sampling)
fixedpars <- fixedparsdata$value
names(fixedpars) <- fixedparsdata$parname

fixed_samples <- list(c(fixedpars, sigma_fixed))

# -----------------------------------------------------------------------------
# Generate initial condition library
# -----------------------------------------------------------------------------
create_lhs_starts <- function(n, lb, ub) {
  if (n <= 0) {
    return(list())
  }
  U <- lhs::randomLHS(n, length(lb))
  log_lb <- log10(lb)
  log_ub <- log10(ub)
  starts <- vector("list", n)
  for (i in seq_len(n)) {
    sample_vec <- 10^(log_lb + U[i, ] * (log_ub - log_lb))
    names(sample_vec) <- fitparnames
    starts[[i]] <- sample_vec
  }
  starts
}

initial_library <- c(list(par_ini), create_lhs_starts(n_random_inits, lower_bounds, upper_bounds))

# Helper tibble describing initial conditions
init_tbl <- tibble(
  start_id = seq_along(initial_library),
  origin = c("baseline", rep("lhs", length(initial_library) - 1))
)

# -----------------------------------------------------------------------------
# Helper to evaluate one fit
# -----------------------------------------------------------------------------
run_one_fit <- function(start_id, solver_name, fixedpars_i, print_level) {
  x0 <- initial_library[[start_id]]
  names(x0) <- fitparnames

  opts <- list(
    algorithm = solver_name,
    maxeval = max_steps,
    maxtime = max_time_seconds,
    print_level = print_level,
    ftol_rel = ftol_rel
  )

  start_time <- Sys.time()
  fit_res <- try(
    nloptr::nloptr(
      x0 = x0,
      eval_f = fit_model_function,
      lb = lower_bounds,
      ub = upper_bounds,
      opts = opts,
      fitdata = fitdata,
      Y0 = Y0,
      tfinal = tfinal,
      dt = dt,
      fitparnames = fitparnames,
      fixedpars = fixedpars_i,
      doses = unique(fitdata$Dose),
      scenarios = scenarios,
      solvertype = solvertype,
      tols = tols
    ),
    silent = TRUE
  )
  runtime <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  if (inherits(fit_res, "try-error")) {
    return(list(
      summary = tibble(
        start_id = start_id,
        solver = solver_name,
        objective = Inf,
        status = NA_integer_,
        iterations = NA_real_,
        runtime_sec = runtime,
        message = as.character(fit_res)
      ),
      bestfit = NULL,
      init = x0
    ))
  }

  solution <- fit_res$solution
  names(solution) <- fitparnames

  bestfit <- list(
    solution = solution,
    objective = fit_res$objective,
    status = fit_res$status,
    iterations = fit_res$iterations %||% NA_real_,
    message = fit_res$message,
    fitpars = setNames(solution, fitparnames),
    fitparnames = fitparnames,
    fixedpars = fixedpars_i,
    Y0 = Y0,
    fitdata = fitdata,
    parlabels = parlabels,
    algorithm = solver_name,
    start_id = start_id,
    parstring = paste0("c(", paste(format(solution, digits = 15), collapse = ", "), ")")
  )

  tibble_row <- tibble(
    start_id = start_id,
    solver = solver_name,
    objective = fit_res$objective,
    status = fit_res$status,
    iterations = fit_res$iterations %||% NA_real_,
    runtime_sec = runtime,
    message = fit_res$message
  )

  list(
    summary = tibble_row,
    bestfit = bestfit,
    init = x0
  )
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# -----------------------------------------------------------------------------
# Main loop over fixed-parameter samples
# -----------------------------------------------------------------------------
summary_list <- list()
bestfit_list <- vector("list", length(fixed_samples))
all_candidates <- vector("list", length(fixed_samples))

for (sample_idx in seq_along(fixed_samples)) {
  fixedpars_i <- fixed_samples[[sample_idx]]

  task_grid <- tidyr::expand_grid(
    start_id = seq_along(initial_library),
    solver = solver_candidates
  )

if (nrow(task_grid) == 1) {
    future::plan(sequential)
    print_level <- 1
  } else {
    workers <- min(max_workers, future::availableCores(), nrow(task_grid))
    future::plan(multisession, workers = workers)
    print_level <- print_level_default
  }

  candidate_results <- future_lapply(
    seq_len(nrow(task_grid)),
    function(idx) {
      task <- task_grid[idx, ]
      run_one_fit(task$start_id, task$solver, fixedpars_i, print_level)
    },
    future.seed = TRUE
  )

  future::plan(sequential)

  candidate_tbl <- bind_rows(purrr::map(candidate_results, "summary")) %>%
    left_join(init_tbl, by = "start_id") %>%
    arrange(objective)

  best_idx <- which.min(candidate_tbl$objective)
  bestfit_entry <- candidate_results[[best_idx]]$bestfit

  if (is.null(bestfit_entry)) {
    stop("All optimisation attempts failed; no finite objective value obtained.")
  }

  bestfit_list[[sample_idx]] <- bestfit_entry
  all_candidates[[sample_idx]] <- list(
    summary = candidate_tbl,
    details = candidate_results
  )
  summary_list[[sample_idx]] <- candidate_tbl %>% mutate(sample = sample_idx)
}

summary_tbl <- bind_rows(summary_list)

# -----------------------------------------------------------------------------
# Persist results
# -----------------------------------------------------------------------------
output_dir <- here::here("results", "output")
fig_dir <- here::here("results", "figures")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

readr::write_csv(summary_tbl, file.path(output_dir, "multistart-summary.csv"))
saveRDS(all_candidates, file.path(output_dir, "multistart-candidates.Rds"))
saveRDS(bestfit_list, file.path(output_dir, "bestfit.Rds"))

# -----------------------------------------------------------------------------
# Generate a best-fit figure for visual inspection
# -----------------------------------------------------------------------------
make_bestfit_plot <- function(bestfit, filename) {
  scenario_map <- c(
    `0` = "NoTreatment",
    `10` = "PanCytoVir10mg",
    `100` = "PanCytoVir100mg"
  )

  # Simulate trajectories for the figure doses (baseline schedule only)
  sim_list <- simulate_dose_predictions(
    bestfit,
    ts_doses = figure_doses,
    solvertype = solvertype,
    tols = tols,
    dt = dt,
    tfinal = tfinal
  )

  ts_df <- sim_list$timeseries_df %>%
    filter(Schedule == "s1", Dose %in% figure_doses) %>%
    mutate(
      Scenario = factor(scenario_map[as.character(Dose)], levels = scenario_map)
    )

  plt <- plot_timeseries(
    data = bestfit$fitdata,
    modelfit = ts_df,
    tmax = tfinal,
    dose_levels = figure_doses,
    dose_levels_labels = figure_dose_labels,
    x_jitter = 0.3
  )

  png(
    filename,
    width = 8,
    height = 5,
    units = "in",
    res = 300
  )
  plot(plt)
  dev.off()
}

if (length(bestfit_list) >= 1 && !is.null(bestfit_list[[1]])) {
  make_bestfit_plot(
    bestfit = bestfit_list[[1]],
    filename = file.path(fig_dir, "bestfit-multistart.png")
  )

  message(
    "Finished multi-start fitting. Best objective: ",
    signif(bestfit_list[[1]]$objective, 6)
  )
} else {
  warning("No successful fits were produced; skipping figure generation.")
}
