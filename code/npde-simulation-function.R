###############################################################################
# npde-simulation-function.R
#
# This file contains helper functions for performing Normalized Prediction
# Distribution Error (NPDE) analyses for the QSP model implemented in this
# project.  Every function is written with extensive inline documentation so
# that readers who are new to R and to NPDE concepts can follow along.
#
# The functions defined here do not load R packages themselves; instead, the
# calling scripts should load the packages they need (dplyr, purrr, tibble,
# readr, etc.).  This mirrors the structure used in the rest of the code base
# where reusable functions remain lightweight and dependency-free.
###############################################################################

#' Load the processed experimental data that were used for model fitting.
#'
#' @param data_path Character string with the path to the processed data CSV.
#'   By default we read the canonical data set stored in
#'   `data/processed-data/processeddata.csv`.
#'
#' @return A tibble with one row per measurement.  Additional helper columns
#'   are added to make downstream processing easier:
#'   * Scenario is stored as an ordered factor (No treatment, 10 mg/kg, 100 mg/kg).
#'   * Quantity identifies the outcome (LogVirusLoad, IL6, or Weight).
#'   * Dose contains the dose in mg/kg for each scenario.
#'   * Day is duplicated as `xvals` to match the notation used elsewhere.
#'
load_observed_dataset <- function(data_path) {
  data <- readr::read_csv(
    file = data_path,
    show_col_types = FALSE
  )

  # Ensure that the scenario ordering matches the order used during model fit.
  scenario_levels <- c("NoTreatment", "PanCytoVir10mg", "PanCytoVir100mg")
  quantity_levels <- c("LogVirusLoad", "IL6", "Weight")

  data <- data %>%
    dplyr::mutate(
      Scenario = factor(Scenario, levels = scenario_levels, ordered = TRUE),
      Quantity = factor(Quantity, levels = quantity_levels, ordered = TRUE),
      # Extract the numeric dose (in mg/kg) from the scenario name.  For the
      # "NoTreatment" scenario the dose is set to zero.
      Dose = dplyr::case_when(
        Scenario == "NoTreatment" ~ 0,
        TRUE ~ readr::parse_number(as.character(Scenario))
      ),
      xvals = Day
    )

  return(data)
}

#' Create a look-up table with dosing information for each scenario.
#'
#' @param observed_data Tibble returned by `load_observed_dataset()`.
#'
#' @return Tibble with two columns: Scenario (character) and dose_mgkg (numeric).
#'
create_scenario_table <- function(observed_data) {
  scenario_df <- observed_data %>%
    dplyr::distinct(Scenario, Dose) %>%
    dplyr::mutate(
      Scenario = as.character(Scenario),
      dose_mgkg = as.numeric(Dose)
    ) %>%
    dplyr::select(Scenario, dose_mgkg) %>%
    # Arrange for reproducibility of downstream operations.
    dplyr::arrange(Scenario)

  return(scenario_df)
}

#' Estimate the measurement variability for each outcome variable.
#'
#' NPDE computations require a description of the residual (measurement) noise.
#' We approximate this by examining the variability among replicate
#' measurements in the experimental data.  For each Quantity we compute the
#' standard deviation across replicates within each (Scenario, Day) group, and
#' then combine those standard deviations using a simple pooled estimate.
#'
#' @param observed_data Tibble returned by `load_observed_dataset()`.
#'
#' @return A tibble with columns `Quantity` and `sd` (one row per quantity).
#'   The `sd` values are guaranteed to be positive; if we cannot estimate an SD
#'   from the data (e.g., only one replicate), we fall back to the overall
#'   standard deviation across all observations, or to a small default value.
#'
estimate_measurement_sd <- function(observed_data) {
  # Compute the replicate-to-replicate standard deviation for each time and
  # scenario.  `stats::sd()` returns NA when only one replicate is available;
  # we retain those NAs for now and handle them later.
  replicate_sd <- observed_data %>%
    dplyr::group_by(Quantity, Scenario, Day) %>%
    dplyr::summarise(
      replicate_sd = stats::sd(Value),
      .groups = "drop"
    )

  # Overall variability for each quantity (used as a fallback when replicate
  # level standard deviations are unavailable).
  overall_sd <- observed_data %>%
    dplyr::group_by(Quantity) %>%
    dplyr::summarise(
      overall_sd = stats::sd(Value),
      .groups = "drop"
    )

  # Combine all available replicate-level SDs into a single pooled estimate for
  # each quantity.  We square the SDs before averaging so that we average
  # variances (which add), and then take the square root to return to SD units.
  pooled_sd <- replicate_sd %>%
    dplyr::filter(!is.na(replicate_sd)) %>%
    dplyr::group_by(Quantity) %>%
    dplyr::summarise(
      pooled_sd = sqrt(mean(replicate_sd^2)),
      .groups = "drop"
    )

  measurement_sd <- overall_sd %>%
    dplyr::left_join(pooled_sd, by = "Quantity") %>%
    dplyr::mutate(
      sd = dplyr::case_when(
        !is.na(pooled_sd) & pooled_sd > 0 ~ pooled_sd,
        !is.na(overall_sd) & overall_sd > 0 ~ overall_sd,
        TRUE ~ 0.1 # conservative default if all else fails
      )
    ) %>%
    dplyr::select(Quantity, sd)

  return(measurement_sd)
}

#' Simulate posterior predictive distributions for each observation time.
#'
#' @param bestfit_list List of posterior samples (each element is a list that
#'   contains the parameter estimates, fixed parameters, and initial conditions).
#' @param scenario_table Output of `create_scenario_table()` describing the
#'   dose for each scenario in mg/kg.
#' @param times Numeric vector of observation times (in days) that we want the
#'   simulator to evaluate.
#' @param tfinal Final simulation time.  Defaults to the maximum of `times` to
#'   avoid unnecessary computation.
#' @param dt Time step used to construct the time grid for which solutions are
#'   returned.  The default of 0.01 days matches the existing workflows.
#' @param txstart,txend,txinterval Dosing schedule parameters.  These default to
#'   the values used during model calibration (start on day 1, dose every 12
#'   hours, stop after day 4).
#'
#' @return A tibble with columns Sample (posterior sample index), Scenario,
#'   Day, Quantity, and Prediction.  Each row corresponds to a single model
#'   prediction (after applying the correct transformation for the quantity).
#'
simulate_posterior_predictions <- function(
  bestfit_list,
  scenario_table,
  times,
  tfinal = max(times),
  dt = 0.01,
  txstart = 1,
  txend = 4,
  txinterval = 0.5
) {
  # Guard against accidental missing inputs that would lead to cryptic errors.
  if (length(bestfit_list) == 0) {
    stop("The list of best-fit samples is empty.  NPDE simulations require at least one sample.")
  }

  if (length(times) == 0) {
    stop("The vector of observation times is empty.  Provide at least one time point.")
  }

  # The simulation end time must be at least as large as the final observation.
  tfinal <- max(tfinal, max(times))

  predictions <- purrr::imap_dfr(
    bestfit_list,
    function(bestfit, sample_id) {
      # Extract the dynamic model inputs from the posterior sample.
      params <- bestfit$solution
      names(params) <- bestfit$fitparnames
      fixedpars <- bestfit$fixedpars
      Y0 <- bestfit$Y0

      scenario_predictions <- purrr::pmap_dfr(
        scenario_table,
        function(Scenario, dose_mgkg) {
          # Convert the mg/kg dose into the "Ad0" parameter by dividing by the
          # 50 g mouse body weight scaling used throughout the project.
          dose_amount <- dose_mgkg / 50

          # Collect all inputs for the simulator function.  Wrapping everything
          # into a named vector keeps the call to `simulate_model()` concise and
          # mirrors the pattern used in other scripts.
          sim_inputs <- c(
            Y0,
            params,
            fixedpars,
            Ad0 = dose_amount,
            txstart = txstart,
            txinterval = txinterval,
            txend = txend,
            tstart = 0,
            tfinal = tfinal,
            dt = dt
          )

          # Run the deterministic simulator.  Errors are allowed to propagate so
          # that the calling environment can diagnose issues explicitly.
          ode_output <- do.call(simulate_model, as.list(sim_inputs))
          ode_df <- tibble::as_tibble(ode_output)

          # Keep only the observation times of interest.  The simulator returns
          # results on a dense grid so we filter by exact matches.  This mirrors
          # how existing scripts work with the ODE output.
          keep_rows <- ode_df %>%
            dplyr::filter(time %in% times)

          if (nrow(keep_rows) == 0) {
            stop("No simulation rows matched the requested observation times.  Check the `times` vector and time step.")
          }

          # Convert state variables to the observable quantities.  Virus load is
          # log-transformed (with a floor at 1 virion), whereas the innate
          # response (F) and symptoms (S) remain on the original scale.
          keep_rows <- keep_rows %>%
            dplyr::mutate(
              Scenario = Scenario,
              Sample = as.integer(sample_id),
              Day = time,
              LogVirusLoad = log10(pmax(1, V)),
              IL6 = F,
              Weight = S
            ) %>%
            dplyr::select(Sample, Scenario, Day, LogVirusLoad, IL6, Weight) %>%
            tidyr::pivot_longer(
              cols = c(LogVirusLoad, IL6, Weight),
              names_to = "Quantity",
              values_to = "Prediction"
            )

          return(keep_rows)
        }
      )

      return(scenario_predictions)
    }
  )

  return(predictions)
}

#' Compute NPDE values for every observation in the data set.
#'
#' @param observed_data Tibble of observed measurements.
#' @param prediction_draws Tibble returned by `simulate_posterior_predictions()`.
#' @param measurement_sd Tibble returned by `estimate_measurement_sd()`.
#'
#' @return Tibble identical to `observed_data` but with additional columns that
#'   describe the predictive distribution and the resulting NPDE value.
#'
compute_npde_table <- function(observed_data, prediction_draws, measurement_sd) {
  # Prepare a named vector so we can easily look up the SD for each quantity
  # while iterating over individual observations.
  measurement_sd_vec <- stats::setNames(measurement_sd$sd, measurement_sd$Quantity)

  # Collect all posterior predictions that correspond to each observation
  # (Scenario, Day, Quantity).  Storing them as list-columns keeps the original
  # values available for descriptive statistics and NPDE calculations.
  prediction_list <- prediction_draws %>%
    dplyr::group_by(Scenario, Day, Quantity) %>%
    dplyr::summarise(
      predictions = list(Prediction),
      .groups = "drop"
    )

  augmented <- observed_data %>%
    dplyr::left_join(prediction_list, by = c("Scenario", "Day", "Quantity"))

  npde_table <- augmented %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      # Retrieve the vector of posterior predictive means.
      prediction_values = list(if (is.null(predictions)) numeric(0) else unlist(predictions)),
      n_predictions = length(prediction_values[[1]]),
      # Compute descriptive statistics of the predictive distribution.  When no
      # predictions are available we return NA to highlight the issue.
      prediction_mean = if (n_predictions > 0) mean(prediction_values[[1]]) else NA_real_,
      prediction_median = if (n_predictions > 0) stats::median(prediction_values[[1]]) else NA_real_,
      prediction_sd = if (n_predictions > 1) stats::sd(prediction_values[[1]]) else NA_real_,
      prediction_lower = if (n_predictions > 0) stats::quantile(prediction_values[[1]], probs = 0.025) else NA_real_,
      prediction_upper = if (n_predictions > 0) stats::quantile(prediction_values[[1]], probs = 0.975) else NA_real_,
      # Measurement SD lookup with a conservative minimum to avoid numerical
      # issues in subsequent probability calculations.
      measurement_sigma = {
        sigma <- measurement_sd_vec[[as.character(Quantity)]]
        if (is.null(sigma) || is.na(sigma) || sigma <= 0) sigma <- 0.1
        sigma
      },
      # Evaluate the predictive cumulative distribution function (CDF) at the
      # observed value by averaging the Normal CDF across posterior samples.
      # Clipping keeps the probability away from the exact 0/1 boundaries so
      # that `qnorm()` remains stable.
      predictive_cdf = if (n_predictions > 0) {
        prob_vals <- stats::pnorm(Value, mean = prediction_values[[1]], sd = measurement_sigma)
        prob_mean <- mean(prob_vals)
        prob_mean <- min(max(prob_mean, .Machine$double.eps), 1 - .Machine$double.eps)
        prob_mean
      } else {
        NA_real_
      },
      NPDE = if (!is.na(predictive_cdf)) stats::qnorm(predictive_cdf) else NA_real_,
      percentile_rank = if (n_predictions > 0) mean(prediction_values[[1]] <= Value) else NA_real_,
      deviation = Value - prediction_mean
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-predictions, -prediction_values, -n_predictions)

  return(npde_table)
}

#' Summarise NPDE diagnostics across each quantity.
#'
#' @param npde_table Output of `compute_npde_table()`.
#'
#' @return Tibble with one row per quantity summarising the mean, standard
#'   deviation, median, Shapiro-Wilk p-value, and the fraction of NPDE values
#'   outside the ±1.96 interval.
#'
create_npde_summary <- function(npde_table) {
  summary_tbl <- npde_table %>%
    dplyr::filter(!is.na(NPDE)) %>%
    dplyr::group_by(Quantity) %>%
    dplyr::summarise(
      mean_npde = mean(NPDE),
      sd_npde = stats::sd(NPDE),
      median_npde = stats::median(NPDE),
      proportion_outside_1_96 = mean(abs(NPDE) > 1.96),
      shapiro_p = if (dplyr::n() >= 3) stats::shapiro.test(NPDE)$p.value else NA_real_,
      .groups = "drop"
    )

  return(summary_tbl)
}

#' Orchestrate the full NPDE analysis pipeline.
#'
#' @param bestfit_list List of posterior samples.
#' @param data_path Path to the processed experimental data (CSV file).
#' @param tfinal Optional simulation end time.
#' @param dt Optional simulation time step.
#'
#' @return A named list that gathers together every intermediate and final
#'   result so that downstream scripts (analysis or plotting) can reuse them.
#'
run_full_npde_analysis <- function(
  bestfit_list,
  data_path,
  tfinal = NULL,
  dt = 0.01
) {
  observed_data <- load_observed_dataset(data_path)
  scenario_table <- create_scenario_table(observed_data)
  observation_times <- sort(unique(observed_data$Day))

  if (is.null(tfinal)) {
    tfinal <- max(observation_times)
  }

  measurement_sd <- estimate_measurement_sd(observed_data)
  prediction_draws <- simulate_posterior_predictions(
    bestfit_list = bestfit_list,
    scenario_table = scenario_table,
    times = observation_times,
    tfinal = tfinal,
    dt = dt
  )
  npde_table <- compute_npde_table(observed_data, prediction_draws, measurement_sd)
  summary_tbl <- create_npde_summary(npde_table)

  result <- list(
    observed_data = observed_data,
    scenario_table = scenario_table,
    measurement_sd = measurement_sd,
    prediction_draws = prediction_draws,
    npde_table = npde_table,
    summary = summary_tbl
  )

  return(result)
}

###############################################################################
# End of file
###############################################################################
