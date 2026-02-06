# -----------------------------------------------------------------------------
# objective-components-function.R
# -----------------------------------------------------------------------------
# PURPOSE
#   Centralized computation of objective function components (NLL, SSR, weights)
#   and standardized residuals. All parts of the codebase should use this file
#   to avoid inconsistent weighting or sigma handling.
# -----------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(here)

# Centralized virus transform helpers.
source(here::here("code", "analysis-code", "functions", "virus-transform-function.R"))

#' Build a long-format prediction table from simulator output.
#'
#' @param sim_df Simulator output with columns time, V, F, S and a scenario column.
#' @param scenario_col Name of the scenario column in sim_df.
#' @param time_col Name of the time column in sim_df.
#' @return Long data.frame with columns: Scenario, Day, Quantity, Predicted.
#'
build_prediction_long <- function(
  sim_df,
  scenario_col = "Scenario",
  time_col = "time"
) {
  if (!all(c(time_col, scenario_col, "V", "F", "S") %in% names(sim_df))) {
    stop("sim_df must contain time, Scenario, V, F, and S columns.")
  }

  virus_col <- virus_quantity_name

  sim_df %>%
    mutate(
      !!virus_col := transform_virus(.data[["V"]]),
      IL6 = .data[["F"]],
      WeightLossPerc = .data[["S"]]
    ) %>%
    select(
      Day = all_of(time_col),
      Scenario = all_of(scenario_col),
      all_of(virus_col),
      IL6,
      WeightLossPerc
    ) %>%
    pivot_longer(
      cols = c(all_of(virus_col), IL6, WeightLossPerc),
      names_to = "Quantity",
      values_to = "Predicted"
    )
}

#' Extract sigma parameters from fitted + fixed parameter vectors.
#'
#' @param params Named vector of fitted parameters (may include sigma_*).
#' @param fixedpars Named vector of fixed parameters (may include sigma_*).
#' @return Named numeric vector of sigma parameters.
#'
collect_sigma_pool <- function(params, fixedpars) {
  vec <- c(params, fixedpars)
  vec <- vec[!is.na(names(vec))]
  vec <- vec[grepl("^sigma_(add|prop)_", names(vec))]

  # De-duplicate by name (keep first occurrence, i.e., fitted over fixed).
  if (length(vec) && any(duplicated(names(vec)))) {
    vec <- vec[!duplicated(names(vec))]
  }

  vec
}

#' Compute objective components and residuals.
#'
#' @param fitdata Standardized fitdata with columns Scenario, Quantity, Day, Value.
#' @param pred_long Long predictions with Scenario, Quantity, Day, Predicted.
#' @param sigma_pool Named vector of sigma parameters (additive + proportional).
#' @param min_variance Minimum variance to avoid log(0) or divide-by-zero.
#' @param time_round Decimal places to round Day/time for safe joins.
#' @param objective_only If TRUE, return only the scalar objective (skip
#'   residuals and breakdown tables) to reduce memory use during fitting.
#' @return List with:
#'   - objective: scalar (sum of weighted NLL blocks)
#'   - breakdown_nll: data.frame with per-block NLL and weights
#'   - breakdown_ssr: data.frame with per-block SSR and weights
#'   - residuals: data.frame with residuals, variance, and weighted residuals
#'
compute_objective_components <- function(
  fitdata,
  pred_long,
  sigma_pool,
  min_variance = 1e-12,
  time_round = 8,
  objective_only = FALSE
) {
  # Validate required columns.
  required_fit <- c("Scenario", "Quantity", "Day", "Value")
  required_pred <- c("Scenario", "Quantity", "Day", "Predicted")
  missing_fit <- setdiff(required_fit, names(fitdata))
  missing_pred <- setdiff(required_pred, names(pred_long))
  if (length(missing_fit)) {
    stop("fitdata is missing required columns: ", paste(missing_fit, collapse = ", "))
  }
  if (length(missing_pred)) {
    stop("pred_long is missing required columns: ", paste(missing_pred, collapse = ", "))
  }

  # Establish stable factor levels.
  scenario_levels <- if (is.factor(fitdata$Scenario)) {
    levels(fitdata$Scenario)
  } else {
    sort(unique(fitdata$Scenario))
  }
  quantity_levels <- if (is.factor(fitdata$Quantity)) {
    levels(fitdata$Quantity)
  } else {
    sort(unique(fitdata$Quantity))
  }

  # Ensure matching types for joins.
  fitdata <- fitdata %>%
    mutate(
      Scenario = factor(Scenario, levels = scenario_levels),
      Quantity = factor(Quantity, levels = quantity_levels)
    )
  pred_long <- pred_long %>%
    mutate(
      Scenario = factor(Scenario, levels = scenario_levels),
      Quantity = factor(Quantity, levels = quantity_levels)
    )

  # Mitigate floating point mismatch between observation days and solver output.
  if (is.numeric(fitdata$Day)) {
    fitdata <- fitdata %>% mutate(Day = round(Day, time_round))
  }
  if (is.numeric(pred_long$Day)) {
    pred_long <- pred_long %>% mutate(Day = round(Day, time_round))
  }

  data_joined <- fitdata %>%
    select(Scenario, Quantity, Day, Value) %>%
    inner_join(pred_long, by = c("Scenario", "Quantity", "Day"))

  if (!nrow(data_joined)) {
    return(list(
      objective = Inf,
      breakdown_nll = NULL,
      breakdown_ssr = NULL,
      residuals = NULL
    ))
  }

  if (nrow(data_joined) < nrow(fitdata)) {
    return(list(
      objective = Inf,
      breakdown_nll = NULL,
      breakdown_ssr = NULL,
      residuals = if (objective_only) NULL else data_joined
    ))
  }

  if (isTRUE(objective_only)) {
    data_joined <- data_joined %>%
      select(Scenario, Quantity, Value, Predicted)
  }

  sigma_pool <- sigma_pool[!is.na(names(sigma_pool))]
  if (length(sigma_pool) && any(duplicated(names(sigma_pool)))) {
    sigma_pool <- sigma_pool[!duplicated(names(sigma_pool))]
  }

  sigma_tbl <- tibble(name = names(sigma_pool), value = as.numeric(sigma_pool)) %>%
    filter(grepl("^sigma_(add|prop)_", name))

  if (!nrow(sigma_tbl)) {
    sigma_tbl <- tibble(
      Quantity = factor(quantity_levels, levels = quantity_levels),
      add = 0,
      prop = 0
    )
  } else {
    sigma_tbl <- sigma_tbl %>%
      tidyr::extract(
        name,
        into = c("type", "Quantity"),
        regex = "^sigma_(add|prop)_(.+)$"
      ) %>%
      select(Quantity, type, value) %>%
      pivot_wider(names_from = type, values_from = value)

    if (!"add" %in% names(sigma_tbl)) sigma_tbl$add <- 0
    if (!"prop" %in% names(sigma_tbl)) sigma_tbl$prop <- 0

    sigma_tbl <- sigma_tbl %>%
      mutate(Quantity = factor(Quantity, levels = quantity_levels))
  }

  weight_tbl <- fitdata %>%
    count(Quantity, Scenario, name = "n") %>%
    mutate(weight = 1 / n)

  resid_df <- data_joined %>%
    left_join(sigma_tbl, by = "Quantity") %>%
    left_join(weight_tbl, by = c("Quantity", "Scenario"))

  if (isTRUE(objective_only)) {
    resid_df <- resid_df %>%
      mutate(
        add = ifelse(is.na(add), 0, add),
        prop = ifelse(is.na(prop), 0, prop),
        Variance = pmax(add^2 + (prop * Predicted)^2, min_variance),
        Residual = Value - Predicted,
        NLLPoint = 0.5 * (log(Variance) + (Residual^2) / Variance)
      )
  } else {
    resid_df <- resid_df %>%
      mutate(
        add = ifelse(is.na(add), 0, add),
        prop = ifelse(is.na(prop), 0, prop),
        Variance = pmax(add^2 + (prop * Predicted)^2, min_variance),
        Residual = Value - Predicted,
        StdResid = Residual / sqrt(Variance),
        NLLPoint = 0.5 * (log(Variance) + (Residual^2) / Variance),
        NLLWeightedResid = sqrt(weight) * StdResid
      )
  }

  if (any(!is.finite(resid_df$NLLPoint))) {
    return(list(
      objective = Inf,
      breakdown_nll = NULL,
      breakdown_ssr = NULL,
      residuals = if (objective_only) NULL else resid_df
    ))
  }

  if (isTRUE(objective_only)) {
    return(list(
      objective = sum(resid_df$NLLPoint * resid_df$weight),
      breakdown_nll = NULL,
      breakdown_ssr = NULL,
      residuals = NULL
    ))
  }

  breakdown_nll <- resid_df %>%
    group_by(Quantity, Scenario, weight) %>%
    summarise(
      n = n(),
      nll = sum(NLLPoint),
      weighted_nll = weight * sum(NLLPoint),
      .groups = "drop"
    )

  breakdown_ssr <- resid_df %>%
    group_by(Quantity, Scenario, weight) %>%
    summarise(
      n = n(),
      ssr = sum(Residual^2),
      weighted_ssr = weight * sum(Residual^2),
      .groups = "drop"
    )

  list(
    objective = sum(breakdown_nll$weighted_nll),
    breakdown_nll = breakdown_nll,
    breakdown_ssr = breakdown_ssr,
    residuals = resid_df
  )
}
