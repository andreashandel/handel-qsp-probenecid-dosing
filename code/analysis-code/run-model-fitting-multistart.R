# R script to fit model1 or model2 with a two-stage, multi-start workflow
#
# -----------------------------------------------------------------------------
# Goal of this script
# -----------------------------------------------------------------------------
# This script provides a single entry point to fit either model1 or model2
# efficiently. It preserves the current model definitions and the existing
# fit_model_function logic, but improves parameter exploration by:
#   1) Sampling many starting points (global exploration) using Latin Hypercube
#      sampling in parameter space.
#   2) Quickly scoring those starts with a coarse ODE setup (fast screening).
#   3) Refining the best candidates with a local optimizer and tighter ODE
#      settings (accurate final fit).
#
# Fixed parameters are loaded once from file and are NOT sampled in this script.
#
# The script keeps a similar flexibility to the existing scripts by allowing
# ad-hoc parameter fixing using `user_fixed_params`.
#
# The code is intentionally verbose and heavily documented for clarity.
#
# -----------------------------------------------------------------------------
# Required packages
# -----------------------------------------------------------------------------
# We use the same packages already referenced in the existing scripts so the
# environment requirements remain familiar.
#
# - here: for robust project-relative file paths
# - dplyr: data manipulation
# - nloptr: local optimization
# - deSolve: ODE solving (used inside the simulator functions)
# - lhs: Latin Hypercube sampling
# - future + future.apply: parallel execution across multiple starts

library(here)
library(dplyr)
library(nloptr)
library(deSolve)
library(lhs)
library(future)
library(future.apply)

# -----------------------------------------------------------------------------
# Source simulator and fit functions (do NOT modify these functions here)
# -----------------------------------------------------------------------------
# The simulator functions define the ODE systems. The fit functions define the
# objective function to be minimized. We keep these as-is to preserve the model.

source(here::here("code", "analysis-code", "model1-simulator-function.R"))
source(here::here("code", "analysis-code", "model2-simulator-function.R"))
source(here::here("code", "analysis-code", "model1-fit-function.R"))
source(here::here("code", "analysis-code", "model2-fit-function.R"))
source(here::here("code", "plotting-code", "stage1-trajectory-plot-function.R"))

# -----------------------------------------------------------------------------
# User settings: choose model and fit options
# -----------------------------------------------------------------------------
# Pick which model to fit. Valid values are "model1" or "model2".
model_choice <- "model1"

# Number of random starts for the global exploration stage.
# Increase this for a more thorough search (but longer runtime).
# Typical values: 20-100 depending on compute budget.
n_starts <- 100

# Stage 1 screening settings.
# The screening can optionally run a *short* local optimization for each start
# to get a better initial objective value. Set stage1_maxeval to a small number
# (e.g., 5-20) to improve screening quality while keeping it fast.
# If stage1_maxeval <= 1, the objective is evaluated once with no optimization.
stage1_algorithm <- "NLOPT_LN_NELDERMEAD"
stage1_maxeval <- 1

# Stage 1 plotting settings.
# Set to TRUE to generate a single multi-panel plot that overlays trajectories
# from all Stage 1 candidates after screening.
# NOTE: This can be slow or produce a very large figure if there are many starts.
# Candidates with objective values above stage1_plot_obj_max are excluded.
stage1_plot_models <- FALSE
stage1_plot_max <- Inf
stage1_plot_file <- here::here("results", "figures", "stage1-candidate-trajectories.png")
stage1_plot_width <- 12
stage1_plot_height <- 16
stage1_plot_dpi <- 300
stage1_plot_dt <- 0.05
stage1_plot_tfinal <- 7
stage1_plot_alpha <- 0.2
stage1_plot_obj_max <- 5000
stage1_plot_max_F <- 10
stage1_plot_max_S <- 50

# Stage 1 candidate filtering for refinement.
# Any candidate with objective > stage1_refine_obj_max will be excluded BEFORE
# applying refine_selection_mode. Use Inf to disable filtering.
stage1_refine_obj_max <- 1e4

# Fraction of prior best fits (from previous run) to include as starting points.
# 0 = include none, 1 = include all, values in-between keep the best fraction.
prior_fit_fraction <- 0.25


# Number of top candidates (best objective values) to refine locally.
# These are selected from the global stage.
n_refine <- 20

# If TRUE, print progress messages during Stage 2 refinement.
stage2_verbose <- TRUE


# How to choose candidates for refinement beyond the single best one.
# Options:
#   "best"   = use the lowest objective values from Stage 1
#   "random" = sample remaining candidates at random
#   "mixed"  = take a mix of best + random (controlled by mix_best_fraction)
#refine_selection_mode <- "best" # "best", "random", or "mixed"
refine_selection_mode <- "mixed" # "best", "random", or "mixed"

# For mixed selection: fraction of remaining slots taken from the best values.
# The rest are chosen at random from the remaining candidates.
mix_best_fraction <- 0.2


# Stage 2 optimizer settings.
# Provide one or more NLOPT algorithm names. Each top candidate is refined with
# EACH algorithm in this vector.
# Examples: "NLOPT_LN_BOBYQA", "NLOPT_LN_NELDERMEAD", "NLOPT_LN_COBYLA"
#stage2_algorithms <- c("NLOPT_LN_BOBYQA", "NLOPT_LN_COBYLA", "NLOPT_LN_SBPLX")
stage2_algorithms <- c("NLOPT_LN_BOBYQA")
#stage2_algorithms <- c("NLOPT_LN_COBYQA")
#stage2_algorithms <- c("NLOPT_LN_SBPLX")

# Parallel workers to use. If NULL, uses available cores.
# Set to 1 to force sequential execution.
#n_workers <- n_refine * length(stage2_algorithms)
n_workers <- NULL

# Maximum number of optimizer iterations for each Stage 2 run.
stage2_maxeval <- 2000

# Fit in log-space for positive parameters? This improves optimizer stability.
# This does NOT change the model; it only changes how we search the space.
use_log_space <- TRUE

# NOTE: The fit_function handles log-space transforms internally (mirrors the
# run-model1-fitting workflow). We pass logfit = 1 when use_log_space is TRUE,
# and keep parameters in log space when calling fit_function.
logfit_in_fit_function <- ifelse(use_log_space, 1, 0)

# Optional: hold some parameters fixed (ad hoc), similar to the existing scripts.
# Provide named vector of parameter values to fix.
# Example: user_fixed_params <- c(Emax_F = 1)
user_fixed_params <- c(Emax_F = 1)

# Random seed control for LHS sampling.
# Use seed_mode = "fixed" for reproducible results, or "time" for a new seed
# each run based on system time.
seed_mode <- "time" # "fixed" or "time"
seed_value <- 1234

if (seed_mode == "fixed") {
  set.seed(seed_value)
} else if (seed_mode == "time") {
  set.seed(as.integer(Sys.time()))
} else {
  stop("seed_mode must be either 'fixed' or 'time'.")
}

# Basic input validation for key settings.
if (!is.numeric(stage1_maxeval) || stage1_maxeval < 1) {
  stop("stage1_maxeval must be a number >= 1.")
}
if (refine_selection_mode == "mixed" &&
    (mix_best_fraction <= 0 || mix_best_fraction >= 1)) {
  stop("mix_best_fraction must be between 0 and 1 when refine_selection_mode = 'mixed'.")
}
if (!length(stage2_algorithms)) {
  stop("stage2_algorithms must contain at least one optimizer name.")
}
if (!is.numeric(stage1_plot_dt) || stage1_plot_dt <= 0) {
  stop("stage1_plot_dt must be a number > 0.")
}
if (!is.numeric(stage1_plot_tfinal) || stage1_plot_tfinal <= 0) {
  stop("stage1_plot_tfinal must be a number > 0.")
}
if (!is.numeric(stage1_plot_obj_max)) {
  stop("stage1_plot_obj_max must be numeric (use Inf for no filtering).")
}
if (!is.numeric(prior_fit_fraction) ||
    prior_fit_fraction < 0 ||
    prior_fit_fraction > 1) {
  stop("prior_fit_fraction must be between 0 and 1.")
}
if (!is.numeric(stage1_refine_obj_max)) {
  stop("stage1_refine_obj_max must be numeric (use Inf for no filtering).")
}
if (!is.numeric(stage1_plot_max_F) || stage1_plot_max_F <= 0) {
  stop("stage1_plot_max_F must be a number > 0.")
}
if (!is.numeric(stage1_plot_max_S) || stage1_plot_max_S <= 0) {
  stop("stage1_plot_max_S must be a number > 0.")
}

# -----------------------------------------------------------------------------
# Helper: load and prepare data shared by both models
# -----------------------------------------------------------------------------
# Load processed data once. We then perform small cleanups to match model
# expectations: Scenario and Quantity factors, Dose calculation, and xvals.

fitdata_path <- here::here("data", "processed-data", "processeddata.csv")
fitdata <- read.csv(fitdata_path, stringsAsFactors = FALSE)

fitdata$Scenario <- factor(
  fitdata$Scenario,
  levels = c("NoTreatment", "PanCytoVir10mg", "PanCytoVir100mg")
)
fitdata$Quantity <- factor(
  fitdata$Quantity,
  levels = c("LogVirusLoad", "IL6", "WeightLossPerc")
)
fitdata$Dose <- c(0, 10, 100)[as.numeric(fitdata$Scenario)]
fitdata$xvals <- fitdata$Day

scenarios <- levels(fitdata$Scenario)
doses <- sort(unique(fitdata$Dose))

# -----------------------------------------------------------------------------
# Helper: build model-specific configuration
# -----------------------------------------------------------------------------
# This function returns all information needed for fitting a given model:
# - initial conditions (Y0)
# - initial parameter values
# - bounds
# - fixed parameter file
# - simulator function

build_model_config <- function(model_choice) {
  if (model_choice == "model1") {
    # Initial conditions for model1 (no E compartment)
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

    # Parameter initial values and bounds copied from run-model1-fitting.R
    par_ini_full <- c(
      b = 1e-8,
      k = 1e-5,
      p = 1e4,
      kF = 0.1,
      cV = 100,
      gF = 1,
      hV = 1e3,
      Fmax = 2,
      hF = 1,
      gS = 10,
      cS = 1,
      Emax_F = 1,
      C50_F = 1e-5,
      C50_V = 1e-8
    )

    lb <- c(
      b = 1e-12,
      k = 1e-10,
      p = 1,
      kF = 1e-1,
      cV = 0.1,
      gF = 1e-3,
      hV = 1e-2,
      Fmax = 0.1,
      hF = 1e-5,
      gS = 1e-3,
      cS = 1e-3,
      Emax_F = 1e-3,
      C50_F = 1e-10,
      C50_V = 1e-10
    )

    ub <- c(
      b = 1e-5,
      k = 1,
      p = 1e7,
      kF = 1e2,
      cV = 1e5,
      gF = 1e3,
      hV = 1e5,
      Fmax = 1e3,
      hF = 1e3,
      gS = 1e3,
      cS = 1e3,
      Emax_F = 1,
      C50_F = 1e2,
      C50_V = 1e2
    )

    fixedpars_file <- here::here("data", "processed-data", "model1-fixed-parameters.csv")
    simulatorname <- model1_simulator
    fit_function <- fit_model1_function

    return(list(
      Y0 = Y0,
      par_ini_full = par_ini_full,
      lb = lb,
      ub = ub,
      fixedpars_file = fixedpars_file,
      simulatorname = simulatorname,
      fit_function = fit_function
    ))
  }

  if (model_choice == "model2") {
    # Initial conditions for model2 (includes E compartment)
    Y0 <- c(
      Ad = 0,
      Ac = 0,
      At = 0,
      U = 1e7,
      E = 0,
      I = 0,
      V = 1,
      F = 0,
      A = 0,
      S = 0
    )

    # Parameter initial values and bounds copied from run-model2-fitting.R
    par_ini_full <- c(
      b = 1e-8,
      k = 1e-5,
      p = 1e4,
      kF = 0.1,
      cV = 100,
      gF = 1,
      hV = 1e3,
      Fmax = 2,
      hF = 1,
      gS = 10,
      cS = 1,
      Emax_F = 1,
      C50_F = 1e-5,
      C50_V = 1e-8
    )

    lb <- c(
      b = 1e-12,
      k = 1e-10,
      p = 1,
      kF = 1e-1,
      cV = 0.1,
      gF = 1e-3,
      hV = 1e-2,
      Fmax = 0.1,
      hF = 1e-5,
      gS = 1e-3,
      cS = 1e-3,
      Emax_F = 1e-3,
      C50_F = 1e-10,
      C50_V = 1e-10
    )

    ub <- c(
      b = 1e-5,
      k = 1,
      p = 1e7,
      kF = 1e2,
      cV = 1e5,
      gF = 1e3,
      hV = 1e5,
      Fmax = 1e3,
      hF = 1e3,
      gS = 1e3,
      cS = 1e3,
      Emax_F = 1,
      C50_F = 1e2,
      C50_V = 1e2
    )

    fixedpars_file <- here::here("data", "processed-data", "model2-fixed-parameters.csv")
    simulatorname <- model2_simulator
    fit_function <- fit_model_function

    return(list(
      Y0 = Y0,
      par_ini_full = par_ini_full,
      lb = lb,
      ub = ub,
      fixedpars_file = fixedpars_file,
      simulatorname = simulatorname,
      fit_function = fit_function
    ))
  }

  stop("Unknown model_choice. Use 'model1' or 'model2'.")
}

# -----------------------------------------------------------------------------
# Load model-specific configuration
# -----------------------------------------------------------------------------
config <- build_model_config(model_choice)

Y0 <- config$Y0
par_ini_full <- config$par_ini_full
lb <- config$lb
ub <- config$ub
simulatorname <- config$simulatorname
fit_function <- config$fit_function

# -----------------------------------------------------------------------------
# Load fixed parameters (model-specific)
# -----------------------------------------------------------------------------
fixedparsdata <- read.csv(config$fixedpars_file, stringsAsFactors = FALSE)
fixedpars <- fixedparsdata[, 3]
names(fixedpars) <- fixedparsdata[, 1]

# -----------------------------------------------------------------------------
# Sigma setup (error model parameters)
# -----------------------------------------------------------------------------
# The objective function supports additive and proportional errors for each
# quantity. Here we compute empirical values and keep all sigmas FIXED.

var_by_qty <- fitdata %>%
  group_by(Quantity) %>%
  summarize(v = var(Value, na.rm = TRUE), .groups = "drop")

sigma_all <- c(
  sigma_add_LogVirusLoad = sqrt(as.numeric(var_by_qty[1, 2])),
  sigma_prop_LogVirusLoad = 0.0,
  sigma_add_IL6 = sqrt(as.numeric(var_by_qty[2, 2])),
  sigma_prop_IL6 = 0.0,
  sigma_add_WeightLossPerc = sqrt(as.numeric(var_by_qty[3, 2])),
  sigma_prop_WeightLossPerc = 0.0
)

# Keep all sigmas fixed at these empirical values.
sigma_fixed <- sigma_all

# -----------------------------------------------------------------------------
# Optional user-fixed parameters (ad hoc)
# -----------------------------------------------------------------------------
# This block mirrors the behavior in the existing scripts:
# - Remove user-fixed parameters from the fit vector
# - Add them later to fixed parameters

if (length(user_fixed_params)) {
  missing_names <- setdiff(names(user_fixed_params), names(par_ini_full))
  if (length(missing_names) > 0) {
    stop("user_fixed_params contains unknown parameter names: ",
         paste(missing_names, collapse = ", "))
  }

  par_ini_full <- par_ini_full[setdiff(names(par_ini_full), names(user_fixed_params))]
  lb <- lb[setdiff(names(lb), names(user_fixed_params))]
  ub <- ub[setdiff(names(ub), names(user_fixed_params))]
}

fitparnames <- names(par_ini_full)

# Add user-fixed parameters to the fixed parameter vector
if (length(user_fixed_params)) {
  fixedpars <- c(fixedpars, user_fixed_params)
}

# Add fixed sigmas to fixed parameter vector
fixedpars <- c(fixedpars, sigma_fixed)

# -----------------------------------------------------------------------------
# Solver settings
# -----------------------------------------------------------------------------
# Use a single set of solver settings for both stages.
# These are the more accurate settings.

solver_settings <- list(
  solvertype = "vode",
  tols = 1e-9,
  tfinal = 7,
  dt = 0.01
)

# -----------------------------------------------------------------------------
# Objective wrapper
# -----------------------------------------------------------------------------
# This function evaluates the objective for a given parameter vector in either
# fast or fine mode. It handles log-space transforms and passes all required
# inputs to fit_model_function.

objective_wrapper <- function(par_vec, solver_settings) {
  # Ensure names are attached to parameter vector
  params <- par_vec
  names(params) <- fitparnames

  # Evaluate objective
  fit_function(
    params = unname(params),
    fitdata = fitdata,
    Y0 = Y0,
    tfinal = solver_settings$tfinal,
    dt = solver_settings$dt,
    fitparnames = fitparnames,
    fixedpars = fixedpars,
    doses = doses,
    scenarios = scenarios,
    solvertype = solver_settings$solvertype,
    tols = solver_settings$tols,
    simulatorname = simulatorname,
    logfit = logfit_in_fit_function
  )
}

# Safe wrapper for the objective so we can capture and report errors.
safe_objective <- function(par_vec, solver_settings) {
  tryCatch(
    list(
      value = objective_wrapper(par_vec, solver_settings),
      error = NULL
    ),
    error = function(e) {
      list(
        value = Inf,
        error = conditionMessage(e)
      )
    }
  )
}

# Stage 1 screening for a single candidate.
# If stage1_maxeval > 1, we run a short optimization to get a better objective
# estimate; otherwise, we evaluate the objective once.
stage1_screen_one <- function(x0) {
  if (stage1_maxeval <= 1) {
    res <- safe_objective(x0, solver_settings)
    res$par <- x0
    return(res)
  }

  tryCatch(
    {
      result <- nloptr::nloptr(
        x0 = as.numeric(x0),
        eval_f = function(x) objective_wrapper(x, solver_settings),
        lb = par_lb,
        ub = par_ub,
        opts = list(
          algorithm = stage1_algorithm,
          maxeval = stage1_maxeval,
          ftol_rel = 1e-6,
          print_level = 0
        )
      )
      list(value = result$objective, error = NULL, par = result$solution)
    },
    error = function(e) {
      list(value = Inf, error = conditionMessage(e), par = NULL)
    }
  )
}

# -----------------------------------------------------------------------------
# Global exploration: generate starting points
# -----------------------------------------------------------------------------
# We create a candidate list in three steps:
#   1) Load previous best fits (if available) and use them as starting points.
#   2) Generate additional candidates using Latin Hypercube sampling (LHS).
#   3) Always include the baseline initial guess as a candidate.
#
# These candidates are then screened quickly (Stage 1) and the best subset is
# refined in Stage 2.

# Work in log space if requested
par_lb <- lb
par_ub <- ub
par_ini_raw <- par_ini_full
par_ini <- par_ini_full

if (use_log_space) {
  par_lb <- log(par_lb)
  par_ub <- log(par_ub)
  par_ini <- log(par_ini)
}

# Helper to extract a start vector from a saved best-fit object.
extract_start_from_bestfit <- function(bestfit, fitparnames, use_log_space, par_ini_raw) {
  if (is.null(bestfit$fitpars)) {
    return(NULL)
  }
  if (is.null(names(bestfit$fitpars))) {
    names(bestfit$fitpars) <- bestfit$fitparnames
  }
  overlap <- intersect(fitparnames, names(bestfit$fitpars))
  if (!length(overlap)) {
    return(NULL)
  }
  par_vec <- par_ini_raw
  par_vec[overlap] <- bestfit$fitpars[overlap]
  if (use_log_space) {
    par_vec <- log(par_vec)
  }
  setNames(as.numeric(par_vec), fitparnames)
}

# Load best fits from previous runs, if available.
previous_fit_file <- here::here(
  "results",
  "output",
  paste0(model_choice, "-bestfit-multistart.Rds")
)

previous_starts <- list()
if (file.exists(previous_fit_file)) {
  previous_fits <- readRDS(previous_fit_file)
  if (is.list(previous_fits)) {
    objectives <- vapply(previous_fits, function(x) {
      if (!is.null(x$objective)) x$objective else Inf
    }, numeric(1))

    order_idx <- order(objectives)
    n_keep <- ceiling(length(previous_fits) * prior_fit_fraction)
    n_keep <- min(n_keep, length(previous_fits))
    if (n_keep > 0) {
      previous_fits <- previous_fits[order_idx[seq_len(n_keep)]]
    } else {
      previous_fits <- list()
    }

    previous_starts <- lapply(previous_fits, extract_start_from_bestfit,
                              fitparnames = fitparnames,
                              use_log_space = use_log_space,
                              par_ini_raw = par_ini_raw)
    previous_starts <- Filter(Negate(is.null), previous_starts)
  }
}

# Build LHS samples in [0,1], then scale to [lb, ub]
if (n_starts > 0) {
  U <- lhs::randomLHS(n_starts, length(par_ini))
  M <- sweep(U, 2, (par_ub - par_lb), "*")
  M <- sweep(M, 2, par_lb, "+")
  colnames(M) <- fitparnames

  start_list <- lapply(seq_len(n_starts), function(i) {
    setNames(as.numeric(M[i, ]), fitparnames)
  })
} else {
  start_list <- list(setNames(as.numeric(par_ini), fitparnames))
}

# Combine previous fits, baseline, and LHS candidates.
start_list <- c(
  previous_starts,
  list(setNames(as.numeric(par_ini), fitparnames)),
  start_list
)

# Remove duplicate starting points (exact matches).
start_df <- do.call(rbind, start_list)
unique_idx <- !duplicated(start_df)
start_list <- start_list[unique_idx]

# Quick sanity check: make sure the baseline evaluation runs at least once.
baseline_check <- safe_objective(setNames(as.numeric(par_ini), fitparnames), solver_settings)
if (is.infinite(baseline_check$value)) {
  stop(
    "Baseline objective evaluation failed. Error: ",
    baseline_check$error
  )
}

# -----------------------------------------------------------------------------
# Parallel setup
# -----------------------------------------------------------------------------
if (is.null(n_workers)) {
  n_workers <- max(1, future::availableCores() - 1)
}

if (n_workers > 1) {
  future::plan(multisession, workers = n_workers)
} else {
  future::plan(sequential)
}

# -----------------------------------------------------------------------------
# Stage 1: Objective screening for all starts
# -----------------------------------------------------------------------------
# We score all starting points to pick the most promising candidates. If
# stage1_maxeval > 1, we run a short local optimization for each start.

message("Stage 1: screening ", length(start_list), " starting points...")

screen_results <- future_lapply(
  start_list,
  function(x0) stage1_screen_one(x0),
  future.seed = TRUE
)

obj_fast <- vapply(screen_results, function(x) x$value, numeric(1))
err_fast <- vapply(screen_results, function(x) {
  if (is.null(x$error)) "" else x$error
}, character(1))
par_stage1 <- lapply(screen_results, function(x) x$par)

# Attach candidate indices to objective values for labeling.
names(obj_fast) <- as.character(seq_along(obj_fast))

if (all(is.infinite(obj_fast))) {
  message("All Stage 1 objectives are Inf. Showing the most common errors:")
  err_table <- sort(table(err_fast), decreasing = TRUE)
  print(err_table[seq_len(min(5, length(err_table)))])
  stop("Stage 1 failed for all candidates. See errors above.")
}

# Identify candidates for refinement.
eligible_idx <- which(is.finite(obj_fast) & obj_fast <= stage1_refine_obj_max)
if (!length(eligible_idx)) {
  stop("No Stage 1 candidates remain after stage1_refine_obj_max filtering.")
}

order_idx <- eligible_idx[order(obj_fast[eligible_idx])]
best_idx <- order_idx[1]
remaining_idx <- setdiff(order_idx, best_idx)

if (n_refine <= 1 || length(remaining_idx) == 0) {
  keep_idx <- best_idx
} else {
  slots_remaining <- min(n_refine - 1, length(remaining_idx))
  if (refine_selection_mode == "best") {
    keep_idx <- c(best_idx, remaining_idx[seq_len(slots_remaining)])
  } else if (refine_selection_mode == "random") {
    random_idx <- sample(remaining_idx, size = slots_remaining, replace = FALSE)
    keep_idx <- c(best_idx, random_idx)
  } else if (refine_selection_mode == "mixed") {
    n_best <- max(1, floor(slots_remaining * mix_best_fraction))
    n_best <- min(n_best, slots_remaining)
    n_rand <- slots_remaining - n_best
    best_part <- remaining_idx[seq_len(n_best)]
    remaining_pool <- setdiff(remaining_idx, best_part)
    if (n_rand > 0 && length(remaining_pool) > 0) {
      rand_part <- sample(remaining_pool, size = n_rand, replace = FALSE)
    } else {
      rand_part <- integer(0)
    }
    keep_idx <- c(best_idx, best_part, rand_part)
  } else {
    stop("refine_selection_mode must be 'best', 'random', or 'mixed'.")
  }
}

start_list_refine <- start_list[keep_idx]
obj_fast_refine <- obj_fast[keep_idx]

message("Stage 1 done. Objective values chosen for refinement: ")
print(obj_fast_refine)

# -----------------------------------------------------------------------------
# Stage 1: Plot trajectories for all screened candidates (optional)
# -----------------------------------------------------------------------------
simulate_candidate <- function(param_vec, candidate_id) {
  if (is.null(param_vec)) {
    return(NULL)
  }

  params <- param_vec
  names(params) <- fitparnames
  if (use_log_space) {
    params <- exp(params)
  }

  params_ode <- params[!grepl("^sigma_", names(params))]
  fixedpars_ode <- fixedpars[!grepl("^sigma_", names(fixedpars))]

  sim_one <- function(ad0, scenario_label) {
    allpars <- c(
      as.list(Y0),
      as.list(params_ode),
      as.list(fixedpars_ode),
      list(
        Ad0 = ad0,
        txstart = 1,
        txinterval = 0.5,
        txend = 3.9,
        tstart = 0,
        tfinal = stage1_plot_tfinal,
        dt = stage1_plot_dt,
        solvertype = solver_settings$solvertype,
        tols = solver_settings$tols
      )
    )

    odeout <- tryCatch(do.call(simulatorname, allpars), error = identity)
    if (inherits(odeout, "error")) {
      return(NULL)
    }

    ode_df <- as.data.frame(odeout)

    # Keep only the required state variables for fast plotting.
    keep_cols <- c("time", "U", "I", "V", "F", "A", "S")
    keep_cols <- keep_cols[keep_cols %in% names(ode_df)]
    ode_df <- ode_df[, keep_cols, drop = FALSE]
    ode_df$Dose <- ad0
    ode_df$Candidate <- candidate_id

    ode_df
  }

  sim_list <- lapply(seq_along(doses), function(i) {
    sim_one(doses[i], scenarios[i])
  })

  sim_list <- Filter(Negate(is.null), sim_list)
  if (!length(sim_list)) {
    return(NULL)
  }

  dplyr::bind_rows(sim_list)
}

if (stage1_plot_models) {
  candidate_ids <- which(is.finite(obj_fast) & obj_fast <= stage1_plot_obj_max)
  if (is.finite(stage1_plot_max)) {
    candidate_ids <- head(candidate_ids, stage1_plot_max)
  }

  if (!length(candidate_ids)) {
    message("No Stage 1 candidates selected for plotting after filtering.")
  } else {
  sim_list <- lapply(candidate_ids, function(i) {
    simulate_candidate(par_stage1[[i]], i)
  })
  sim_list <- Filter(Negate(is.null), sim_list)

  if (length(sim_list)) {
    sim_all <- dplyr::bind_rows(sim_list)

    # Ensure figures folder exists.
    dir.create(here::here("results", "figures"), showWarnings = FALSE, recursive = TRUE)

    stage1_plot <- plot_stage1_trajectories(
      sim_df = sim_all,
      data = fitdata,
      candidate_objectives = obj_fast,
      output_file = stage1_plot_file,
      width = stage1_plot_width,
      height = stage1_plot_height,
      dpi = stage1_plot_dpi,
      show_plot = FALSE,
      alpha_lines = stage1_plot_alpha,
      max_F = stage1_plot_max_F,
      max_S = stage1_plot_max_S
    )

    # Load the saved figure from file and display it.
    if (!requireNamespace("png", quietly = TRUE)) {
      stop("Package 'png' is required to display the saved plot file.")
    }
    img <- png::readPNG(stage1_plot_file)
    if (!interactive()) {
      grDevices::dev.new()
    }
    grid::grid.newpage()
    grid::grid.raster(img)
    grDevices::dev.flush()

    message("Stage 1 plot saved to: ", stage1_plot_file)
    user_choice <- readline("Stage 1 plot displayed. Press [enter] to continue or type 'q' to stop: ")
    if (tolower(trimws(user_choice)) == "q") {
      stop("Execution stopped by user after Stage 1 plot.")
    }
  } else {
    message("No Stage 1 plot produced (all simulations failed).")
  }
  }
}

# for debugging
#start_df <- do.call(rbind, start_list)

#
#browser()

# -----------------------------------------------------------------------------
# Stage 2: Local refinement of top candidates
# -----------------------------------------------------------------------------
# We now run a local optimizer starting from each of the best points.
# Each candidate is refined with EACH algorithm listed in stage2_algorithms.
# This uses the shared solver_settings for all runs.

local_optimize_one <- function(x0, algorithm) {
  # nloptr expects numeric vector; names are stored separately
  x0_vec <- as.numeric(x0)

  result <- nloptr::nloptr(
    x0 = x0_vec,
    eval_f = function(x) objective_wrapper(x, solver_settings),
    lb = par_lb,
    ub = par_ub,
    opts = list(
      algorithm = algorithm,
      maxeval = stage2_maxeval,
      ftol_rel = 1e-10,
      print_level = 0
    )
  )

  # Attach names to solution
  sol <- result$solution
  names(sol) <- fitparnames

  # Convert back from log space if needed
  if (use_log_space) {
    sol <- exp(sol)
  }

  result$fitpars <- sol
  result$fitparnames <- fitparnames
  result$fixedpars <- fixedpars
  result$Y0 <- Y0
  result$fitdata <- fitdata
  result$simulator <- model_choice
  result$algorithm <- algorithm

  return(result)
}

stage2_grid <- expand.grid(
  candidate_index = seq_along(start_list_refine),
  algorithm = stage2_algorithms,
  stringsAsFactors = FALSE
)

message(
  sprintf(
    "Stage 2: refining top candidates with optimizers: %s | candidates: %d | total runs: %d | maxeval per run: %d",
    paste(stage2_algorithms, collapse = ", "),
    length(start_list_refine),
    nrow(stage2_grid),
    stage2_maxeval
  )
)

bestfits <- future_lapply(
  seq_len(nrow(stage2_grid)),
  function(i) {
    candidate_index <- stage2_grid$candidate_index[i]
    algorithm <- stage2_grid$algorithm[i]
    if (stage2_verbose) {
      message(
        sprintf(
          "Refinement %d/%d started at %s (candidate %d, algorithm %s)",
          i,
          nrow(stage2_grid),
          format(Sys.time(), "%H:%M:%S"),
          candidate_index,
          algorithm
        )
      )
    }

    start_obj <- obj_fast_refine[candidate_index]
    result <- local_optimize_one(start_list_refine[[candidate_index]], algorithm)

    if (stage2_verbose) {
      message(
        sprintf(
          "Refinement %d/%d finished at %s (start objective = %.6g, final objective = %.6g, optimizer = %s, steps = %d)",
          i,
          nrow(stage2_grid),
          format(Sys.time(), "%H:%M:%S"),
          start_obj,
          result$objective,
          result$algorithm,
          if (!is.null(result$iterations)) result$iterations else NA_integer_
        )
      )
    }

    result
  },
  future.seed = TRUE
)

# -----------------------------------------------------------------------------
# Select best result and save
# -----------------------------------------------------------------------------
# We pick the fit with the smallest objective value.

best_objectives <- vapply(bestfits, function(x) x$objective, numeric(1))
best_idx <- which.min(best_objectives)
bestfit <- bestfits[[best_idx]]

message("Best objective after refinement: ", bestfit$objective)

# Save results for later plotting or reuse
output_file <- here::here(
  "results",
  "output",
  paste0(model_choice, "-bestfit-multistart.Rds")
)

saveRDS(bestfits, output_file)
message("Saved all refined fits to: ", output_file)

# -----------------------------------------------------------------------------
# Cleanup: reset future plan to sequential
# -----------------------------------------------------------------------------
future::plan(sequential)

# End of script
