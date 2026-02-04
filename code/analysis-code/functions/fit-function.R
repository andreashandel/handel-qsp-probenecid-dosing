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

  # Build a sigma pool that includes both FITTED and FIXED sigma parameters.
  sigma_pool <- c(params[fit_sigmas], fixedpars[fixed_sigmas])

  # Accumulate predictions across all dosing scenarios.
  pred_long_list <- list()

  # If we fit in log-space, transform parameters back to natural scale.
  if (logfit == 1) {
    fitpars_ode <- exp(fitpars_ode)
  }

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

    # Extract observed data for each quantity at this scenario.
    Vvals <- fitdata %>%
      filter(Quantity == "LogVirusLoad", Scenario == scenarios[i])
    Innvals <- fitdata %>%
      filter(Quantity == "IL6", Scenario == scenarios[i])
    Symvals <- fitdata %>%
      filter(Quantity == "WeightLossPerc", Scenario == scenarios[i])

    # Simulator time vector for aligning predictions to observation times.
    tvec <- odeout[, "time"]

    # Model predictions aligned to observed times.
    Vpred = log10(pmax(1, odeout[match(Vvals$xvals, tvec), "V"]))
    Innpred = odeout[match(Innvals$xvals, tvec), "F"]
    Sympred = odeout[match(Symvals$xvals, tvec), "S"]

    pred_long_list[[length(pred_long_list) + 1]] <- data.frame(
      Scenario = scenarios[i],
      Day = c(Vvals$xvals, Innvals$xvals, Symvals$xvals),
      Quantity = c(
        rep("LogVirusLoad", length(Vvals$xvals)),
        rep("IL6", length(Innvals$xvals)),
        rep("WeightLossPerc", length(Symvals$xvals))
      ),
      Predicted = c(Vpred, Innpred, Sympred)
    )
  }

  pred_long <- bind_rows(pred_long_list)
  components <- compute_objective_components(fitdata, pred_long, sigma_pool)

  if (!is.finite(components$objective)) {
    return(Inf)
  }

  # Return total objective value to be minimized by the optimizer.
  return(components$objective)
}
