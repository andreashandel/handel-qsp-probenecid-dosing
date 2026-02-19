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
  logfit,
  obs_times_by_scenario = NULL,
  objective_data = NULL,
  weight_mode = "per_block"
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

  if (is.null(obs_times_by_scenario)) {
    if (!is.null(objective_data)) {
      obs_times_by_scenario <- lapply(
        objective_data$data_by_scenario,
        function(x) sort(unique(c(0, x$time)))
      )
    } else {
      obs_times_by_scenario <- split(fitdata$Day, fitdata$Scenario)
      obs_times_by_scenario <- lapply(obs_times_by_scenario, function(x) {
        sort(unique(c(0, x)))
      })
    }
  }

  if (is.null(objective_data) && !weight_mode %in% c("per_block", "equal_quantity")) {
    stop("weight_mode must be 'per_block' or 'equal_quantity'.")
  }

  if (!is.null(objective_data)) {
    quantity_levels <- objective_data$quantity_levels
    quantity_id_map <- objective_data$quantity_id_map

    sigma_add_by_id <- rep(0, length(quantity_levels))
    sigma_prop_by_id <- rep(0, length(quantity_levels))
    for (i in seq_along(quantity_levels)) {
      qname <- quantity_levels[i]
      add_name <- paste0("sigma_add_", qname)
      prop_name <- paste0("sigma_prop_", qname)
      if (!is.null(names(sigma_pool)) && add_name %in% names(sigma_pool)) {
        sigma_add_by_id[i] <- sigma_pool[[add_name]]
      }
      if (!is.null(names(sigma_pool)) && prop_name %in% names(sigma_pool)) {
        sigma_prop_by_id[i] <- sigma_pool[[prop_name]]
      }
    }

    virus_id <- quantity_id_map[[virus_quantity_name]]
    il6_id <- quantity_id_map[["IL6"]]
    wl_id <- quantity_id_map[["WeightLossPerc"]]
  }

  # Loop over all treatment scenarios (one dose per scenario).
  for (i in seq_along(doses)) {
    scenario_key <- as.character(scenarios[i])
    scenario_data <- if (!is.null(objective_data)) {
      objective_data$data_by_scenario[[scenario_key]]
    } else {
      NULL
    }
    if (!is.null(objective_data) && is.null(scenario_data)) {
      next
    }

    obs_times <- obs_times_by_scenario[[scenario_key]]
    if (is.null(obs_times) || length(obs_times) == 0) {
      obs_times <- if (is.null(scenario_data)) 0 else sort(unique(c(0, scenario_data$time)))
    }

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
        tols = tols,
        times = obs_times
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

    # Treat premature termination as a failed fit attempt.
    # Compare against the last requested output time for this scenario.
    target_time <- max(obs_times)
    if (!nrow(ode_df) || max(ode_df$time, na.rm = TRUE) < (target_time - dt / 2)) {
      return(1e10)
    }

    if (!is.null(objective_data)) {
      time_index <- match(scenario_data$time, ode_df$time)
      if (any(is.na(time_index))) {
        return(1e10)
      }

      n_obs <- length(time_index)
      pred <- numeric(n_obs)
      qid <- scenario_data$quantity_id

      if (!is.null(virus_id) && any(qid == virus_id)) {
        v_pred <- transform_virus(ode_df$V[time_index])
        pred[qid == virus_id] <- v_pred[qid == virus_id]
      }
      if (!is.null(il6_id) && any(qid == il6_id)) {
        f_pred <- ode_df$F[time_index]
        pred[qid == il6_id] <- f_pred[qid == il6_id]
      }
      if (!is.null(wl_id) && any(qid == wl_id)) {
        s_pred <- ode_df$S[time_index]
        pred[qid == wl_id] <- s_pred[qid == wl_id]
      }

      add <- sigma_add_by_id[qid]
      prop <- sigma_prop_by_id[qid]
      variance <- pmax(add^2 + (prop * pred)^2, objective_data$min_variance)
      residual <- scenario_data$value - pred
      nll_point <- 0.5 * (log(variance) + (residual^2) / variance)

      if (any(!is.finite(nll_point))) {
        return(1e10)
      }

      objective_total <- objective_total + sum(scenario_data$weight * nll_point)
    } else {
      ode_df$Scenario <- scenarios[i]
      pred_long <- build_prediction_long(ode_df, scenario_col = "Scenario", time_col = "time")
      fitdata_i <- fitdata[fitdata$Scenario == scenarios[i], , drop = FALSE]
      components <- compute_objective_components(
        fitdata_i,
        pred_long,
        sigma_pool,
        weight_mode = weight_mode,
        objective_only = TRUE
      )

      if (!is.finite(components$objective)) {
        return(1e10)
      }

      objective_total <- objective_total + components$objective
    }
  }

  # Return total objective value to be minimized by the optimizer.
  return(objective_total)
}
