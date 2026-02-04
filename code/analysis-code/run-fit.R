# -----------------------------------------------------------------------------
# run-fit.R
# -----------------------------------------------------------------------------
# PURPOSE
#   Unified entry point for model fitting. This script ALWAYS runs the
#   multi-start workflow first (global screening + local refinement). After
#   that, it can optionally run a single-fit workflow that samples fixed
#   parameters while seeding EVERY sample from the multistart bestfit.
#
# OUTPUTS (results/output/)
#   Multistart stage (always):
#     - <model>-bestfit-multistart.Rds
#   Optional fixed-parameter sampling stage:
#     - <model>-bestfit-sample.Rds
# -----------------------------------------------------------------------------
#
# DETAILED WALKTHROUGH
#   - "User settings": choose model + stage toggles, then configure the relevant
#     settings (multistart + optional sampling).
#   - "Common setup": load data and model config, compute sigma settings, and
#     prepare fixed parameters.
#   - "Multistart stage": generate candidate starts, screen them, then refine,
#     and save the multistart bestfits.
#   - "Sampling stage (optional)": run single fits for fixed-parameter samples
#     using the multistart bestfit as the starting point for every sample.
#   - "Save results": write RDS files and report elapsed time.
# -----------------------------------------------------------------------------

rm(list = ls(all.names = TRUE)) # Start with a clean workspace.

library(here)         # Project-root-relative file paths.
library(nloptr)       # Optimizer used for parameter fitting.
library(lhs)          # Latin Hypercube Sampling for sampling/starts.
library(deSolve)      # ODE solver interface.
library(future)       # Parallel plan setup.
library(future.apply) # Parallel lapply for multi-sample fits.

# Model simulators + fit function
source(here::here("code", "analysis-code", "functions", "model1-simulator-function.R"))
source(here::here("code", "analysis-code", "functions", "model2-simulator-function.R"))
source(here::here("code", "analysis-code", "functions", "fit-function.R"))

# Shared helpers
source(here::here("code", "analysis-code", "functions", "fit-data-function.R"))
source(here::here("code", "analysis-code", "functions", "fixed-parameters-function.R"))
source(here::here("code", "analysis-code", "functions", "sigma-settings-function.R"))
source(here::here("code", "analysis-code", "functions", "model-config-function.R"))
source(here::here("code", "analysis-code", "functions", "fixed-parameter-sampling-function.R"))
source(here::here("code", "analysis-code", "functions", "bestfit-pack-function.R"))
source(here::here("code", "analysis-code", "functions", "fit-single-function.R"))

# -----------------------------------------------------------------------------
# User settings (common to all stages)
# -----------------------------------------------------------------------------
model_choice <- "model2" # "model1" or "model2".

# If FALSE, only the multistart stage is run (no fixed-parameter sampling).
run_sampling_stage <- TRUE

# Parallel workers for both stages.
#   - NULL = let each stage choose a sensible default based on available cores.
#   - Any positive integer will cap workers in both stages.
n_workers <- NULL

# Out-of-bounds handling for initial parameter values (both stages).
# Options:
#   "stop"  = abort run if any initial values fall outside bounds
#   "clamp" = replace out-of-bounds values with nearest bound and continue
oob_action <- "clamp"

# Sigma settings (keep fixed by default)
sigma_to_fit <- character(0)

# Optional ad-hoc fixed parameters (removed from fit vector)
user_fixed_params <- c(Emax_F = 1)

# ODE solver settings (used by both stages)
solver_settings <- list(
  solvertype = "vode",
  tols = 1e-8,
  tfinal = 7,
  dt = 0.05
)

# -----------------------------------------------------------------------------
# User settings (multistart stage)
# -----------------------------------------------------------------------------
# NOTE: We only save bestfit-multistart and bestfit-sample files.

# Number of random starts for global exploration.
n_starts <- 200

# Stage 1 screening settings (fast local search or single evaluation).
#   - stage1_maxeval = 1 means "evaluate objective once at the start point".
#   - stage1_maxeval > 1 runs a very short local search to de-noise poor starts.
stage1_algorithm <- "NLOPT_LN_COBYLA"
stage1_maxeval <- 1

# Candidate filtering after Stage 1 screening.
#   - We compute an objective value for each candidate.
#   - Any candidate with objective > stage1_refine_obj_max is discarded.
#   - This is a "hard cutoff" to avoid spending refinement time on clearly
#     implausible parameter sets.
stage1_refine_obj_max <- 1e4

# Prior multistart best fits to include as starting points.
#   - If a previous multistart file exists, we can reuse its best solutions
#     as additional candidates for the current run.
#   - prior_fit_fraction defines the fraction of previous fits to keep
#     (ordered from best to worst objective).
prior_fit_fraction <- 1

# Number of top candidates to refine.
n_refine <- 30

# Refinement selection mode for Stage 2.
#   - "best": choose the top n_refine candidates strictly by Stage 1 objective.
#   - "random": choose n_refine candidates uniformly at random.
#   - "mixed": always include the best candidate, then select a mix of the
#     next-best and random candidates. mix_best_fraction controls how much of
#     the remaining slots are filled by the next-best candidates.
refine_selection_mode <- "mixed" # "best", "random", or "mixed"
mix_best_fraction <- 0.2

# Stage 2 optimizers.
#   - Every candidate selected for Stage 2 is refined with EACH algorithm
#     listed here (cartesian product). This makes results more robust to
#     optimizer-specific quirks but increases runtime proportionally.
#   - Useful Options: "NLOPT_LN_COBYLA", "NLOPT_LN_SBPLX",
#     "NLOPT_LN_NELDERMEAD", "NLOPT_LN_BOBYQA"
#stage2_algorithms <- c("NLOPT_LN_COBYLA", "NLOPT_LN_BOBYQA")
#stage2_algorithms <- c("NLOPT_LN_COBYLA")
#stage2_algorithms <- c("NLOPT_LN_SBPLX")
stage2_algorithms <- c("NLOPT_LN_BOBYQA")

# Stage 2 optimizer settings.
stage2_maxeval <- 500

# Fit in log space for positive parameters (multistart mode).
#   - TRUE means the optimizer sees log(parameters), which ensures positivity
#     and often improves scaling when parameters span orders of magnitude.
use_log_space <- TRUE

# Random seed settings for multistart sampling.
seed_mode <- "time" # "fixed" or "time"
seed_value <- 1234

# -----------------------------------------------------------------------------
# User settings (sampling stage: single-fit mode for fixed-parameter samples)
# -----------------------------------------------------------------------------
# Number of random fixed-parameter samples. 0 = only baseline fixed parameters.
nsamp <- 100

# Force specific fixed parameters to a chosen value across all samples.
fixed_overrides <- c(Emax_V = 1)

# Optimizer settings for the sampling stage.
#   - sample_logfit applies the same log-parameter transformation used in
#     multistart, but only for the sampling stage.
#sample_algorithm <- "NLOPT_LN_BOBYQA"
sample_algorithm <- "NLOPT_LN_COBYLA"
sample_maxeval <- 2000
sample_ftol_rel <- 1e-10
sample_logfit <- TRUE

# If nsamp == 0 (single fit), set sample_print_level > 0 for optimizer output.
sample_print_level <- 1

# -----------------------------------------------------------------------------
# Validate optimizer algorithm names (fail fast with a clear message)
# -----------------------------------------------------------------------------
# This prevents late failures inside nloptr when an invalid algorithm name
# is supplied in user settings.
valid_nloptr_algorithms <- c(
  "NLOPT_GN_DIRECT",
  "NLOPT_GN_DIRECT_L",
  "NLOPT_GN_DIRECT_L_RAND",
  "NLOPT_GN_DIRECT_NOSCAL",
  "NLOPT_GN_DIRECT_L_NOSCAL",
  "NLOPT_GN_DIRECT_L_RAND_NOSCAL",
  "NLOPT_GN_ORIG_DIRECT",
  "NLOPT_GN_ORIG_DIRECT_L",
  "NLOPT_GD_STOGO",
  "NLOPT_GD_STOGO_RAND",
  "NLOPT_LD_SLSQP",
  "NLOPT_LD_LBFGS",
  "NLOPT_LN_PRAXIS",
  "NLOPT_LD_VAR1",
  "NLOPT_LD_VAR2",
  "NLOPT_LD_TNEWTON",
  "NLOPT_LD_TNEWTON_RESTART",
  "NLOPT_LD_TNEWTON_PRECOND",
  "NLOPT_LD_TNEWTON_PRECOND_RESTART",
  "NLOPT_GN_CRS2_LM",
  "NLOPT_GN_MLSL",
  "NLOPT_GD_MLSL",
  "NLOPT_GN_MLSL_LDS",
  "NLOPT_GD_MLSL_LDS",
  "NLOPT_LD_MMA",
  "NLOPT_LD_CCSAQ",
  "NLOPT_LN_COBYLA",
  "NLOPT_LN_NEWUOA",
  "NLOPT_LN_NEWUOA_BOUND",
  "NLOPT_LN_NELDERMEAD",
  "NLOPT_LN_SBPLX",
  "NLOPT_LN_AUGLAG",
  "NLOPT_LD_AUGLAG",
  "NLOPT_LN_AUGLAG_EQ",
  "NLOPT_LD_AUGLAG_EQ",
  "NLOPT_LN_BOBYQA",
  "NLOPT_GN_ISRES",
  "NLOPT_GN_ESCH"
)

validate_algorithms <- function(algos, label) {
  bad <- setdiff(algos, valid_nloptr_algorithms)
  if (length(bad) > 0) {
    stop(
      label,
      " contains invalid NLOPT algorithms: ",
      paste(bad, collapse = ", "),
      ". Valid options include: ",
      paste(valid_nloptr_algorithms, collapse = ", ")
    )
  }
}

validate_algorithms(sample_algorithm, "sample_algorithm")
validate_algorithms(stage1_algorithm, "stage1_algorithm")
validate_algorithms(stage2_algorithms, "stage2_algorithms")

# -----------------------------------------------------------------------------
# Execution timing and run header
# -----------------------------------------------------------------------------
start_time_wall <- Sys.time()
message("Starting run-fit for model: ", model_choice)

# -----------------------------------------------------------------------------
# Common setup (data + model config)
# -----------------------------------------------------------------------------
fitdata_bundle <- load_fit_data()
fitdata <- fitdata_bundle$fitdata
scenarios <- fitdata_bundle$scenarios
doses <- fitdata_bundle$doses

config <- build_model_config(model_choice)

Y0 <- config$Y0
par_ini_full <- config$par_ini_full
lb <- config$lb
ub <- config$ub
parlabels_full <- config$parlabels_full
fixedpars_file <- config$fixedpars_file
simulatorname <- config$simulatorname
fit_function <- config$fit_function

# -----------------------------------------------------------------------------
# Sigma configuration (shared)
# -----------------------------------------------------------------------------
sigma_settings <- compute_sigma_settings(fitdata, sigma_to_fit)

# Append fitted sigma parameters to the fitted parameter vector.
par_ini_full <- c(par_ini_full, sigma_settings$sigma_fit_ini)

lb <- c(lb, rep(1e-6, length(sigma_settings$sigma_fit_ini)))
ub <- c(ub, rep(1e3, length(sigma_settings$sigma_fit_ini)))

names(lb) <- names(par_ini_full)
names(ub) <- names(par_ini_full)

# -----------------------------------------------------------------------------
# Optional user-fixed parameters (shared)
# -----------------------------------------------------------------------------
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

parlabels <- parlabels_full[fitparnames]
if (any(is.na(parlabels))) {
  stop("Parlabel definitions missing entries for some fitted parameters.")
}
if (length(parlabels) != length(par_ini_full)) {
  stop("length of parlabels does not match length of par_ini_full")
}

# -----------------------------------------------------------------------------
# Fixed parameters (shared)
# -----------------------------------------------------------------------------
fixedpars_bundle <- load_fixed_parameters(fixedpars_file)
fixedpars <- fixedpars_bundle$values

if (length(user_fixed_params)) {
  fixedpars <- c(fixedpars, user_fixed_params)
}

fixedpars <- c(fixedpars, sigma_settings$sigma_fixed)

# ---------------------------------------------------------------------------
# Check initial values against bounds (shared)
# ---------------------------------------------------------------------------
handle_oob_initials <- function(par_ini, lb, ub, action_label) {
  oob_low <- par_ini < lb
  oob_high <- par_ini > ub
  if (any(oob_low | oob_high)) {
    bad_names <- names(par_ini)[oob_low | oob_high]
    warning(
      sprintf(
        "%s: initial values out of bounds for: %s",
        action_label,
        paste(bad_names, collapse = ", ")
      )
    )
    if (oob_action == "stop") {
      stop("Stopping due to out-of-bounds initial values (oob_action = 'stop').")
    }
    if (oob_action == "clamp") {
      par_ini <- pmax(pmin(par_ini, ub), lb)
    } else {
      stop("oob_action must be 'stop' or 'clamp'.")
    }
  }
  par_ini
}

# Apply bounds check to baseline initial values before fitting.
par_ini_full <- handle_oob_initials(par_ini_full, lb, ub, "Baseline")

# -----------------------------------------------------------------------------
# Multistart stage (global screening + local refinement)
# -----------------------------------------------------------------------------
{
  if (seed_mode == "fixed") {
    set.seed(seed_value)
  } else if (seed_mode == "time") {
    set.seed(as.integer(Sys.time()))
  } else {
    stop("seed_mode must be either 'fixed' or 'time'.")
  }

  logfit_in_fit_function <- ifelse(use_log_space, 1, 0)

  objective_wrapper <- function(par_vec, solver_settings) {
    params <- par_vec
    names(params) <- fitparnames

    obj <- fit_function(
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
    if (!is.finite(obj)) {
      return(Inf)
    }
    obj
  }

  safe_objective <- function(par_vec, solver_settings) {
    tryCatch(
      list(
        value = objective_wrapper(par_vec, solver_settings),
        error = NULL
      ),
      error = function(e) {
        list(value = Inf, error = conditionMessage(e))
      }
    )
  }

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

  par_lb <- lb
  par_ub <- ub
  par_ini_raw <- par_ini_full
  par_ini <- par_ini_full

  if (use_log_space) {
    par_lb <- log(par_lb)
    par_ub <- log(par_ub)
    par_ini <- log(par_ini)
  }

  # Ensure baseline x0 is within bounds (especially after log transform).
  if (any(par_ini < par_lb | par_ini > par_ub)) {
    if (oob_action == "stop") {
      stop("Stopping due to out-of-bounds initial values (oob_action = 'stop').")
    }
    if (oob_action == "clamp") {
      par_ini <- pmax(pmin(par_ini, par_ub), par_lb)
    } else {
      stop("oob_action must be 'stop' or 'clamp'.")
    }
  }

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
    # Clamp to bounds in natural space to avoid invalid x0.
    par_vec <- pmax(pmin(par_vec, ub), lb)
    if (use_log_space) {
      par_vec <- log(par_vec)
    }
    setNames(as.numeric(par_vec), fitparnames)
  }

  # If a previous multistart run exists, include its best solutions as
  # additional starting points. This improves reproducibility and can reduce
  # the time to rediscover a known good region of parameter space.
  previous_fit_file <- here::here(
    "results",
    "output",
    paste0(model_choice, "-bestfit-multistart.Rds")
  );

  previous_starts <- list();
  if (file.exists(previous_fit_file)) {
    previous_fits <- readRDS(previous_fit_file);
    if (is.list(previous_fits)) {
      objectives <- vapply(previous_fits, function(x) {
        if (!is.null(x$objective)) x$objective else Inf
      }, numeric(1));

      order_idx <- order(objectives);
      n_keep <- ceiling(length(previous_fits) * prior_fit_fraction);
      n_keep <- min(n_keep, length(previous_fits));
      if (n_keep > 0) {
        previous_fits <- previous_fits[order_idx[seq_len(n_keep)]];
      } else {
        previous_fits <- list();
      }

      previous_starts <- lapply(
        previous_fits,
        extract_start_from_bestfit,
        fitparnames = fitparnames,
        use_log_space = use_log_space,
        par_ini_raw = par_ini_raw
      );
      previous_starts <- Filter(Negate(is.null), previous_starts);
    }
  }

  # Create new random starting points (Latin Hypercube in the transformed
  # parameter space). These are exploratory candidates for Stage 1.
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

  # Combine:
  #   1) prior best-fit starts (if any),
  #   2) the baseline start derived from par_ini_full,
  #   3) random LHS starts.
  # The baseline start is always included so there is at least one known
  # sensible candidate even if LHS draws are poor.
  start_list <- c(
    previous_starts,
    list(setNames(as.numeric(par_ini), fitparnames)),
    start_list
  )

  start_df <- do.call(rbind, start_list)
  unique_idx <- !duplicated(start_df)
  start_list <- start_list[unique_idx]

  # Sanity check: the baseline start must yield a finite objective.
  # If it fails, the model/data configuration is inconsistent and we stop.
  baseline_check <- safe_objective(setNames(as.numeric(par_ini), fitparnames), solver_settings)
  if (is.infinite(baseline_check$value)) {
    stop("Baseline objective evaluation failed. Error: ", baseline_check$error)
  }

  workers_multistart <- if (is.null(n_workers)) {
    max(1, future::availableCores() - 1)
  } else {
    n_workers
  }
  workers_multistart <- min(workers_multistart, future::availableCores())

  if (workers_multistart > 1) {
    future::plan(multisession, workers = workers_multistart)
  } else {
    future::plan(sequential)
  }

  message(
    "Stage 1: screening ",
    length(start_list),
    " starting points with up to ",
    workers_multistart,
    " worker(s)..."
  )

  screen_results <- future_lapply(
    start_list,
    function(x0) stage1_screen_one(x0),
    future.seed = TRUE
  )

  obj_fast <- vapply(screen_results, function(x) x$value, numeric(1))
  err_fast <- vapply(screen_results, function(x) {
    if (is.null(x$error)) "" else x$error
  }, character(1))

  if (all(is.infinite(obj_fast))) {
    message("All Stage 1 objectives are Inf. Showing the most common errors:")
    err_table <- sort(table(err_fast), decreasing = TRUE)
    print(err_table[seq_len(min(5, length(err_table)))])
    stop("Stage 1 failed for all candidates. See errors above.")
  }

  # Apply the Stage 1 cutoff rule:
  #   - keep candidates with finite objectives
  #   - and objectives <= stage1_refine_obj_max
  # This prevents Stage 2 from wasting time on clearly bad parameter sets.
  eligible_idx <- which(is.finite(obj_fast) & obj_fast <= stage1_refine_obj_max)
  if (!length(eligible_idx)) {
    stop("No Stage 1 candidates remain after stage1_refine_obj_max filtering.")
  }

  order_idx <- eligible_idx[order(obj_fast[eligible_idx])]
  best_idx <- order_idx[1]
  remaining_idx <- setdiff(order_idx, best_idx)

  # Decide which candidates go to Stage 2.
  #   - The best candidate from Stage 1 is always included.
  #   - The remaining slots are chosen according to refine_selection_mode.
  if (n_refine <= 1 || length(remaining_idx) == 0) {
    keep_idx <- best_idx
  } else {
    slots_remaining <- min(n_refine - 1, length(remaining_idx))
    if (refine_selection_mode == "best") {
      # Take the next-best candidates after the best one.
      keep_idx <- c(best_idx, remaining_idx[seq_len(slots_remaining)])
    } else if (refine_selection_mode == "random") {
      # Sample uniformly from remaining candidates, ignoring rank.
      random_idx <- sample(remaining_idx, size = slots_remaining, replace = FALSE)
      keep_idx <- c(best_idx, random_idx)
    } else if (refine_selection_mode == "mixed") {
      # Mixed mode: take some top-ranked and some random candidates.
      # mix_best_fraction controls the fraction of slots assigned to top-ranked
      # candidates; the remainder are random to preserve diversity.
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

  local_optimize_one <- function(x0, algorithm) {
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

    pack_bestfit(
      fit_result = result,
      fitparnames = fitparnames,
      fixedpars = fixedpars,
      Y0 = Y0,
      fitdata = fitdata,
      parlabels = parlabels,
      algorithm = algorithm,
      logfit = use_log_space
    )
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

      start_obj <- obj_fast_refine[candidate_index]
      run_start <- Sys.time()
      result <- local_optimize_one(start_list_refine[[candidate_index]], algorithm)
      run_end <- Sys.time()
      run_minutes <- as.numeric(difftime(run_end, run_start, units = "mins"))
      steps_taken <- if (!is.null(result$iterations)) result$iterations else NA_integer_

      message(
        sprintf(
          "Refinement %d/%d finished (start objective = %.6g, final objective = %.6g, optimizer = %s, steps = %s, time = %.2f min)",
          i,
          nrow(stage2_grid),
          start_obj,
          result$objective,
          result$algorithm,
          ifelse(is.na(steps_taken), "NA", as.character(steps_taken)),
          run_minutes
        )
      )

      result
    },
    future.seed = TRUE
  )

  best_objectives <- vapply(bestfits, function(x) x$objective, numeric(1))
  best_idx <- which.min(best_objectives)
  bestfit <- bestfits[[best_idx]]
  # Reorder so bestfit is first, followed by all others (preserve structure).
  bestfits <- c(list(bestfit), bestfits[-best_idx])

  message("Best objective after refinement: ", bestfit$objective)

  output_dir <- here::here("results", "output")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  output_file <- here::here(
    "results",
    "output",
    paste0(model_choice, "-bestfit-multistart.Rds")
  )

  saveRDS(bestfits, output_file)
  message("Saved all refined fits to: ", output_file)

  future::plan(sequential)
}

# -----------------------------------------------------------------------------
# Optional sampling stage (single-fit workflow for fixed-parameter samples)
# -----------------------------------------------------------------------------
if (isTRUE(run_sampling_stage)) {
  # Build the list of fixed-parameter sets to evaluate.
  #   - The first element is always the baseline fixed-parameter set.
  #   - The remaining nsamp elements are random perturbations around baseline.
  fixed_samples <- make_fixed_parameter_samples(fixedpars, nsamp, seed = 1234)

  # Optionally override specific fixed parameters across all samples.
  # This is useful for "what-if" analyses where certain values are held fixed.
  fixed_samples <- lapply(fixed_samples, function(x) {
    if (length(fixed_overrides)) {
      x[names(fixed_overrides)] <- fixed_overrides
    }
    x
  })

  sample_bestfit_file <- paste0(model_choice, "-bestfit-sample.Rds")
  sample_bestfit_path <- here::here("results", "output", sample_bestfit_file)

  multistart_path <- here::here(
    "results",
    "output",
    paste0(model_choice, "-bestfit-multistart.Rds")
  )

  # Load the multistart bestfit; this is the starting point for EVERY sample.
  if (!file.exists(multistart_path)) {
    stop(
      "Multistart bestfit file not found at ",
      multistart_path,
      ". Run the multistart stage first."
    )
  }

  multistart_list <- readRDS(multistart_path)
  if (!is.list(multistart_list) || length(multistart_list) == 0) {
    stop("Multistart bestfit file is empty; cannot seed sampling stage.")
  }
  multistart_bestfit <- multistart_list[[1]]
  message("Seeding sampling stage from multistart bestfit.")

  bestfit_all <- vector("list", length(fixed_samples))
  start_time_cpu <- proc.time()

  fit_one_sample <- function(i, print_level = 0) {
    message(sprintf("processing sample %d at %s", i, format(Sys.time(), "%H:%M:%S")))

    fixedpars_i <- fixed_samples[[i]]

    # Choose the starting point for this sample:
    #   - Always use the multistart bestfit so each sample starts
    #     from the same best-fitting parameter vector.
    previous_bestfit <- multistart_bestfit

    run_single_fit(
      fit_function = fit_function,
      par_ini_full = par_ini_full,
      lb = lb,
      ub = ub,
      fitparnames = fitparnames,
      fixedpars = fixedpars_i,
      Y0 = Y0,
      fitdata = fitdata,
      doses = doses,
      scenarios = scenarios,
      simulatorname = simulatorname,
      logfit = sample_logfit,
      algorithm = sample_algorithm,
      maxeval = sample_maxeval,
      ftol_rel = sample_ftol_rel,
      tols = solver_settings$tols,
      solvertype = solver_settings$solvertype,
      tfinal = solver_settings$tfinal,
      dt = solver_settings$dt,
      parlabels = parlabels,
      print_level = print_level,
      oob_action = oob_action,
      previous_bestfit = previous_bestfit
    )
  }

  if (length(fixed_samples) > 1) {
    workers_sampling <- if (is.null(n_workers)) {
      future::availableCores()
    } else {
      n_workers
    }
    workers_sampling <- min(workers_sampling, future::availableCores())
    future::plan(multisession, workers = workers_sampling)
    message("Running sampling stage in parallel with ", workers_sampling, " workers.")

    bestfit_all <- future_lapply(
      seq_along(fixed_samples),
      fit_one_sample,
      print_level = 0,
      future.seed = TRUE
    )

    future::plan(sequential)
  } else {
    message("Single sample detected; running sequentially.")
    bestfit_all[[1]] <- fit_one_sample(1, print_level = sample_print_level)
  }

  end_time_cpu <- proc.time() - start_time_cpu
  runtime_minutes <- end_time_cpu[[3]] / 60
  cat("sampling stage fit took this many minutes:", runtime_minutes, "\n")
  cat("************** \n")
  cat("used algorithm:", sample_algorithm, "\n")
  cat("************** \n")

  new_objectives <- vapply(bestfit_all, function(x) x$objective, numeric(1))
  objective_summary <- data.frame(
    new_objective = new_objectives
  )

  print(objective_summary)

  output_dir <- here::here("results", "output")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  if (file.exists(sample_bestfit_path)) {
    file.copy(
      from = sample_bestfit_path,
      to = here::here("results", "output", paste0("old", sample_bestfit_file)),
      overwrite = TRUE
    )
  }
  saveRDS(bestfit_all, sample_bestfit_path)

  message("Saved sample bestfit list to: ", sample_bestfit_path)
} else {
  message("Sampling stage disabled; only multistart fit was run.")
}

# -----------------------------------------------------------------------------
# Final timing report
# -----------------------------------------------------------------------------
end_time_wall <- Sys.time()
elapsed_minutes <- as.numeric(difftime(end_time_wall, start_time_wall, units = "mins"))
message("run-fit finished for model: ", model_choice)
message(sprintf("Total elapsed time: %.2f minutes", elapsed_minutes))
