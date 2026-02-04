# -----------------------------------------------------------------------------
# sigma-settings-function.R
# -----------------------------------------------------------------------------
# PURPOSE
#   Compute measurement-error (sigma) parameters from the processed data and
#   split them into fitted vs. fixed subsets.
#
# WHY THIS FILE EXISTS
#   Multiple workflows (fitting, diagnostics, Shiny app) rely on the same sigma
#   defaults. Centralizing the calculation ensures consistency.
# -----------------------------------------------------------------------------

library(dplyr) # Grouped variance summaries for each quantity.

#' Compute sigma parameters and split into fitted vs fixed subsets.
#'
#' The variance for each measurement type is estimated empirically from the
#' observed data, then transformed into sigma values. Users can choose to fit
#' some sigma parameters and fix the rest.
#'
#' @param fitdata Processed fitting data returned by load_fit_data().
#' @param sigma_to_fit Character vector of sigma parameter names to fit. Any
#'   sigma not listed here will be treated as fixed.
#' @return A list with:
#'   - sigma_all: named vector with all sigma values
#'   - sigma_fit_ini: subset to be fitted
#'   - sigma_fixed: subset to be fixed
#'
compute_sigma_settings <- function(fitdata, sigma_to_fit = character(0)) {
  # Empirical variance per measurement type; used as sigma defaults.
  var_by_qty <- fitdata %>%
    group_by(Quantity) %>%
    summarize(v = var(Value, na.rm = TRUE), .groups = "drop")

  # Build the full sigma vector (additive + proportional terms).
  sigma_all <- c(
    sigma_add_LogVirusLoad = sqrt(as.numeric(var_by_qty[1, 2])),
    sigma_prop_LogVirusLoad = 0.0,
    sigma_add_IL6 = sqrt(as.numeric(var_by_qty[2, 2])),
    sigma_prop_IL6 = 0.0,
    sigma_add_WeightLossPerc = sqrt(as.numeric(var_by_qty[3, 2])),
    sigma_prop_WeightLossPerc = 0.0
  )

  # Split sigma parameters into fitted vs. fixed subsets.
  sigma_fit_ini <- sigma_all[sigma_to_fit]
  sigma_fixed <- sigma_all[setdiff(names(sigma_all), sigma_to_fit)]

  # Return everything so downstream scripts can choose what they need.
  list(
    sigma_all = sigma_all,
    sigma_fit_ini = sigma_fit_ini,
    sigma_fixed = sigma_fixed
  )
}
