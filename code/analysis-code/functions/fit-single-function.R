# -----------------------------------------------------------------------------
# fit-single-function.R
# -----------------------------------------------------------------------------
# PURPOSE
#   Provide a single, well-documented wrapper around nloptr::nloptr so both the
#   single-fit and multistart pipelines can reuse the same optimization call.
# -----------------------------------------------------------------------------

library(nloptr) # Nonlinear optimization routines.

#' Run a single optimization for one fixed-parameter set.
#'
#' @param fit_function Objective function (model1 or model2).
#' @param par_ini_full Named vector of initial values for fitted parameters.
#' @param lb Named vector of lower bounds.
#' @param ub Named vector of upper bounds.
#' @param fitparnames Character vector of fitted parameter names.
#' @param fixedpars Named vector of fixed parameters (incl. fixed sigmas).
#' @param Y0 Initial condition vector.
#' @param fitdata Standardized fit data.
#' @param doses Numeric doses from fitdata.
#' @param scenarios Scenario names from fitdata.
#' @param simulatorname Simulator function.
#' @param logfit Logical, TRUE if we fit in log space.
#' @param algorithm NLOPT algorithm name.
#' @param maxeval Maximum optimizer iterations.
#' @param ftol_rel Relative tolerance.
#' @param tols ODE solver tolerances.
#' @param solvertype ODE solver type (e.g., "vode" or "lsoda").
#' @param tfinal Final time for the ODE solver.
#' @param dt Time step for the ODE solver.
#' @param obs_times_by_scenario Optional named list of observation times per scenario.
#' @param objective_data Optional precomputed objective data for fast evaluation.
#' @param weight_mode Weighting scheme for observations ("per_block" or "equal_quantity").
#' @param parlabels Parameter labels for output tables.
#' @param print_level nloptr print level (0 = silent, higher = more output).
#' @param oob_action How to handle out-of-bounds initial values ("stop" or "clamp").
#' @param previous_bestfit Optional prior bestfit list to initialize parameters.
#' @return Bestfit object (nloptr result with added metadata).
#'
run_single_fit <- function(
  fit_function,
  par_ini_full,
  lb,
  ub,
  fitparnames,
  fixedpars,
  Y0,
  fitdata,
  doses,
  scenarios,
  simulatorname,
  logfit,
  algorithm,
  maxeval,
  ftol_rel,
  tols,
  solvertype,
  tfinal,
  dt,
  obs_times_by_scenario = NULL,
  objective_data = NULL,
  weight_mode = "per_block",
  parlabels,
  print_level = 0,
  oob_action = "stop",
  previous_bestfit = NULL
) {
  # If a prior best-fit exists, use it as the starting point where possible.
  if (!is.null(previous_bestfit) && !is.null(previous_bestfit$fitpars)) {
    replace_idx <- intersect(names(par_ini_full), names(previous_bestfit$fitpars))
    par_ini_full[replace_idx] <- previous_bestfit$fitpars[replace_idx]
  }

  # Handle out-of-bounds initial values before transforming/logging.
  oob_low <- par_ini_full < lb
  oob_high <- par_ini_full > ub
  if (any(oob_low | oob_high)) {
    bad_names <- names(par_ini_full)[oob_low | oob_high]
    warning(
      sprintf(
        "Initial values out of bounds for: %s",
        paste(bad_names, collapse = ", ")
      )
    )
    if (oob_action == "stop") {
      stop("Stopping due to out-of-bounds initial values (oob_action = 'stop').")
    }
    if (oob_action == "clamp") {
      par_ini_full <- pmax(pmin(par_ini_full, ub), lb)
    } else {
      stop("oob_action must be 'stop' or 'clamp'.")
    }
  }

  # Convert named vectors to numeric vectors for nloptr input.
  par_ini <- as.numeric(par_ini_full)
  lb_used <- lb
  ub_used <- ub

  # Transform the optimization variables to log-space if requested.
  if (isTRUE(logfit)) {
    par_ini <- log(par_ini)
    lb_used <- log(lb_used)
    ub_used <- log(ub_used)
  }

  # Final safeguard: ensure x0 is within bounds (handles log-space rounding).
  if (any(par_ini < lb_used | par_ini > ub_used)) {
    if (oob_action == "stop") {
      stop("Stopping due to out-of-bounds x0 after transforms (oob_action = 'stop').")
    }
    if (oob_action == "clamp") {
      par_ini <- pmax(pmin(par_ini, ub_used), lb_used)
    } else {
      stop("oob_action must be 'stop' or 'clamp'.")
    }
  }

  # Run the optimizer. The fit_function receives all extra arguments listed here.
  fit_result <- nloptr::nloptr(
    x0 = par_ini,
    eval_f = fit_function,
    lb = lb_used,
    ub = ub_used,
    opts = list(
      algorithm = algorithm,
      maxeval = maxeval,
      ftol_rel = ftol_rel,
      print_level = print_level
    ),
    fitdata = fitdata,
    Y0 = Y0,
    tfinal = tfinal,
    dt = dt,
    fitparnames = fitparnames,
    fixedpars = fixedpars,
    doses = doses,
    scenarios = scenarios,
    solvertype = solvertype,
    tols = tols,
    simulatorname = simulatorname,
    logfit = ifelse(logfit, 1, 0),
    obs_times_by_scenario = obs_times_by_scenario,
    objective_data = objective_data,
    weight_mode = weight_mode
  )

  # Attach standardized metadata to the raw nloptr result.
  pack_bestfit(
    fit_result = fit_result,
    fitparnames = fitparnames,
    fixedpars = fixedpars,
    Y0 = Y0,
    fitdata = fitdata,
    parlabels = parlabels,
    algorithm = algorithm,
    logfit = logfit
  )
}
