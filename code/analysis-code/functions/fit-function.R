###################################################################
# fit-function.R
# -----------------------------------------------------------------
# PURPOSE
#   Unified objective function for fitting either model1 or model2.
#
# HOW IT WORKS
#   - The simulator is passed in via `simulatorname`.
#   - The data are compared to predictions on the observation scale.
#   - Measurement error uses additive + proportional sigma parameters.
###################################################################

library(nloptr)  # Optimization interface (used by nloptr::nloptr).
library(dplyr)   # Data manipulation (filter, group_by, summarize).
library(deSolve) # ODE solver utilities (used by simulator functions).
library(here)    # Project-root-relative paths for sourcing helpers.

# Centralized objective component computations (shared across codebase).
source(here::here("code", "analysis-code", "functions", "objective-components-function.R"))

# This function is set up to allow some parameters to be fitted (stored in params)
# and some parameters to be fixed (stored in fixedpars).
fit_model_function <- function(
  params,
  fitdata,
  Y0,
  tfinal,
  dt,
  fitparnames,
  fixedpars,
  doses,
  scenarios,
  solvertype,
  tols,
  simulatorname,
  logfit
) {
  # Guard against malformed inputs early.
  if (any(!is.finite(params)) || any(!is.finite(fixedpars))) {
    return(Inf)
  }

  # nloptr strips names from parameters; restore names so we can index safely.
  names(params) = fitparnames

  # Separate fitted parameters into ODE parameters vs. sigma parameters.
  fit_sigmas <- grepl("^sigma_(add|prop)_", names(params))
  fitpars_ode = params[!fit_sigmas]

  # Separate fixed parameters into ODE parameters vs. sigma parameters.
  fixed_sigmas <- grepl("^sigma_(add|prop)_", names(fixedpars))
  fixedpars_ode = fixedpars[!fixed_sigmas]

  # Accumulate objective contributions across scenarios to avoid holding
  # full simulation outputs in memory at once.
  objective_total <- 0

  # If we fit in log-space, transform parameters back to natural scale.
  if (logfit == 1) {
    fitpars_ode <- exp(fitpars_ode)
    if (any(fit_sigmas)) {
      params[fit_sigmas] <- exp(params[fit_sigmas])
    }
  }

  # Build a sigma pool that includes both FITTED and FIXED sigma parameters.
  sigma_pool <- c(params[fit_sigmas], fixedpars[fixed_sigmas])

  # Loop over all treatment scenarios (one dose per scenario).
  for (i in seq_along(doses)) {
    # Combine parameters together to be sent to the simulator function.
    # Treatment start/end/intervals are the same for all scenarios.
    allpars <- c(
      as.list(Y0),
      as.list(fitpars_ode),
      as.list(fixedpars_ode),
      list(
        Ad0 = doses[i],
        txstart = 1,
        txinterval = 0.5,
        txend = 3.9,
        tstart = 0,
        tfinal = tfinal,
        dt = dt,
        solvertype = solvertype,
        tols = tols
      )
    )

    # Call the simulator; guard against integrator failures.
    odeout <- try(do.call(simulatorname, allpars), silent = TRUE)
    if (inherits(odeout, "try-error")) {
      cat("!!!!!!unresolvable integrator error - triggering early return from optimizer!!!!!!")
      return(1e10)
    }

    ode_df <- as.data.frame(odeout)
    keep_cols <- c("time", "V", "F", "S")
    ode_df <- ode_df[, keep_cols, drop = FALSE]
    ode_df$Scenario <- scenarios[i]

    pred_long <- build_prediction_long(ode_df, scenario_col = "Scenario", time_col = "time")
    fitdata_i <- fitdata[fitdata$Scenario == scenarios[i], , drop = FALSE]
    components <- compute_objective_components(
      fitdata_i,
      pred_long,
      sigma_pool,
      objective_only = TRUE
    )

    if (!is.finite(components$objective)) {
      return(1e10)
    }

    objective_total <- objective_total + components$objective
  }

  # Return total objective value to be minimized by the optimizer.
  return(objective_total)
}
