# -----------------------------------------------------------------------------
# dose-predictions-function.R
# -----------------------------------------------------------------------------
# PURPOSE
#   Simulate model trajectories and dose-response summaries across a range of
#   dosing schedules. This function is shared by model1 and model2.
# -----------------------------------------------------------------------------

library(deSolve) # ODE solver interface.
library(dplyr)   # Data manipulation for summarization.
library(caTools) # trapz() for AUC calculations.
library(here)

# Centralized virus transform helpers.
source(here::here("code", "analysis-code", "functions", "virus-transform-function.R"))

#' Simulate dose-response outcomes and time-series for a fitted model.
#'
#' @param bestfit A single bestfit object (one list element from the fit output).
#' @param ts_doses Numeric vector of doses for which full time-series should be
#'   retained (all other doses only contribute summary AUCs).
#' @param all_doses Numeric vector of ALL doses to simulate (must include ts_doses).
#' @param schedule_defs List of dosing schedules. Each element should be a list
#'   with txstart, txend, txinterval, name, and optional label.
#' @param solvertype ODE solver name (e.g., "lsoda" or "vode").
#' @param tols Solver tolerance.
#' @param dt Time step for output.
#' @param tfinal Final time for simulation.
#' @param simulatorname Simulator function to call.
#' @return A list with summary tables and full trajectories.
#'
simulate_dose_predictions <- function(
  bestfit,
  ts_doses,
  all_doses,
  schedule_defs,
  solvertype,
  tols,
  dt,
  tfinal,
  simulatorname
) {
  if (!all(ts_doses %in% all_doses)) {
    missing_ts <- setdiff(ts_doses, all_doses)
    stop("all_doses must include all ts_doses. Missing: ", paste(missing_ts, collapse = ", "))
  }
  all_doses <- sort(unique(all_doses))

  # ---------------------------------------------------------------------------
  # Inner helper: compute percent reduction relative to baseline in each schedule
  # ---------------------------------------------------------------------------
  compute_percent_reduction <- function(df) {
    # Use Dose == 0 as the baseline for each schedule to avoid reliance on row order.
    df %>%
      group_by(Schedule) %>%
      mutate(
        base_AUCV = AUCV[Dose == 0][1],
        base_AUCF = AUCF[Dose == 0][1],
        base_AUCS = AUCS[Dose == 0][1],
        perc_AUCV = (base_AUCV - AUCV) / base_AUCV * 100,
        perc_AUCF = (base_AUCF - AUCF) / base_AUCF * 100,
        perc_AUCS = (base_AUCS - AUCS) / base_AUCS * 100
      ) %>%
      select(-base_AUCV, -base_AUCF, -base_AUCS)
  }

  # ---------------------------------------------------------------------------
  # Inner helper: simulate one dosing schedule across many doses
  # ---------------------------------------------------------------------------
  simulate_dose_response <- function(
    all_doses,
    txstart,
    txend,
    txinterval,
    schedule_name,
    ts_doses,
    Y0,
    fitpars_ode,
    fixedpars_ode,
    solvertype,
    tols,
    dt,
    tfinal
  ) {
    time_tolerance <- dt / 2
    output_times <- seq(0, tfinal, by = dt)
    # Store the AUC summaries for each dose in this schedule.
    summary_df <- tibble(
      Dose = all_doses,
      AUCV = NA_real_,
      AUCF = NA_real_,
      AUCS = NA_real_,
      Schedule = schedule_name
    )

    # Store full trajectories only for the requested doses.
    ts_store <- list()
    failures <- list()
    warnings <- list()

    for (i in seq_along(all_doses)) {
      # Assemble the full parameter list for the simulator.
      allpars <- c(
        as.list(Y0),
        as.list(fitpars_ode),
        as.list(fixedpars_ode),
        list(
          Ad0 = all_doses[i],
          txstart = txstart,
          txinterval = txinterval,
          txend = txend,
          tstart = 0,
          tfinal = tfinal,
          dt = dt,
          solvertype = solvertype,
          tols = tols,
          times = output_times
        )
      )

      # Run the ODE simulator and capture warnings/errors.
      odeout <- withCallingHandlers(
        tryCatch(
          do.call(simulatorname, allpars),
          error = function(e) e
        ),
        warning = function(w) {
          warnings[[length(warnings) + 1]] <- tibble(
            Schedule = schedule_name,
            Dose = all_doses[i],
            Message = conditionMessage(w)
          )
          invokeRestart("muffleWarning")
        }
      )
      if (inherits(odeout, "error")) {
        failures[[length(failures) + 1]] <- tibble(
          Schedule = schedule_name,
          Dose = all_doses[i],
          Message = conditionMessage(odeout)
        )
        next
      }

      ode_df <- as.data.frame(odeout)
      if (!all(is.finite(ode_df$time)) ||
          any(!is.finite(ode_df$V)) ||
          any(!is.finite(ode_df$F)) ||
          any(!is.finite(ode_df$S))) {
        failures[[length(failures) + 1]] <- tibble(
          Schedule = schedule_name,
          Dose = all_doses[i],
          Message = "Non-finite values in ODE output."
        )
        next
      }
      if (max(ode_df$time, na.rm = TRUE) < (tfinal - time_tolerance)) {
        failures[[length(failures) + 1]] <- tibble(
          Schedule = schedule_name,
          Dose = all_doses[i],
          Message = "ODE output ended before tfinal; skipping incomplete trajectory."
        )
        next
      }

      # Add free drug and PD effect variables for plotting; mirror ODE definitions.
      Ct <- ode_df$At / as.numeric(allpars["Vt"]) # Total drug concentration.
      fu <- as.numeric(allpars["fmax"]) * Ct / (as.numeric(allpars["f50"]) + Ct) # Free fraction.
      Cu <- fu * Ct # Free drug concentration.
      fV <- as.numeric(allpars["Emax_V"]) * Cu / (as.numeric(allpars["C50_V"]) + Cu) # Effect on virus.
      fF <- as.numeric(allpars["Emax_F"]) * Cu / (as.numeric(allpars["C50_F"]) + Cu) # Effect on innate.

      # Append derived variables to the trajectory output.
      ode_df$Cu <- Cu
      ode_df$fV <- fV
      ode_df$fF <- fF

      # AUCs for outcomes of interest.
      summary_df$AUCV[i] <- caTools::trapz(
        ode_df$time,
        transform_virus(ode_df$V)
      )
      summary_df$AUCF[i] <- caTools::trapz(ode_df$time, ode_df$F) # Innate response AUC.
      summary_df$AUCS[i] <- caTools::trapz(ode_df$time, ode_df$S) # Symptom AUC.

      # Add metadata columns for schedule and dose.
      traj_df <- dplyr::mutate(
        ode_df,
        Dose = all_doses[i],
        Schedule = schedule_name
      )

      # Store trajectories only for selected doses to limit output size.
      if (all_doses[i] %in% ts_doses) {
        ts_store[[length(ts_store) + 1]] <- traj_df
      }
    }

    # Return both the summary table and the time-series data.
    list(
      summary = summary_df,
      timeseries = if (length(ts_store)) bind_rows(ts_store) else tibble(),
      failures = if (length(failures)) bind_rows(failures) else tibble(),
      warnings = if (length(warnings)) bind_rows(warnings) else tibble()
    )
  }

  # ---------------------------------------------------------------------------
  # Prepare model parameters for simulation
  # ---------------------------------------------------------------------------
  # Extract fitted and fixed parameters from the bestfit object.
  fitpars <- bestfit$fitpars
  fixedpars <- bestfit$fixedpars
  Y0 <- bestfit$Y0

  # Separate ODE parameters from sigma parameters.
  fit_sigmas <- grepl("^sigma_(add|prop)_", names(fitpars))
  fitpars_ode <- fitpars[!fit_sigmas]

  # Remove sigma parameters from the fixed parameter vector.
  fixed_sigmas <- grepl("^sigma_(add|prop)_", names(fixedpars))
  fixedpars_ode <- fixedpars[!fixed_sigmas]

  # ---------------------------------------------------------------------------
  # Run all schedules
  # ---------------------------------------------------------------------------
  # Run each schedule and collect summary + time-series outputs.
  all_results <- lapply(
    schedule_defs,
    function(s) {
      simulate_dose_response(
        all_doses,
        s$txstart,
        s$txend,
        s$txinterval,
        s$name,
        ts_doses,
        Y0,
        fitpars_ode,
        fixedpars_ode,
        solvertype,
        tols,
        dt,
        tfinal
      )
    }
  )

  # Combine schedule-specific outputs into a single data.frame each.
  all_results_df <- bind_rows(lapply(all_results, `[[`, "summary"))
  timeseries_df <- bind_rows(lapply(all_results, `[[`, "timeseries"))
  failures_df <- bind_rows(lapply(all_results, `[[`, "failures"))
  warnings_df <- bind_rows(lapply(all_results, `[[`, "warnings"))

  if (nrow(failures_df)) {
    warning(
      "Dose-response simulations had failures for ",
      nrow(failures_df),
      " dose-schedule combinations. These entries were skipped."
    )
  }
  if (nrow(warnings_df)) {
    warning(
      "Dose-response simulations emitted ",
      nrow(warnings_df),
      " warnings. See the warnings table in the output object."
    )
  }

  # Compute percent reduction relative to baseline in each schedule.
  label_map <- vapply(
    schedule_defs,
    function(s) if (!is.null(s$label)) s$label else s$name,
    character(1)
  )
  reduction_df <- compute_percent_reduction(all_results_df) %>%
    mutate(Scenario = label_map[Schedule])

  # Return a structured list used by plotting scripts.
  list(
    all_results_df = all_results_df,
    reduction_df = reduction_df,
    timeseries_df = timeseries_df,
    ts_doses = ts_doses,
    failures = failures_df,
    warnings = warnings_df
  )
}
