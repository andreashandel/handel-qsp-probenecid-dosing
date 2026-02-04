# -----------------------------------------------------------------------------
# fit-data-function.R
# -----------------------------------------------------------------------------
# PURPOSE
#   Provide a single, well-documented way to load the processed data and apply
#   the exact factor ordering and column conventions required by the fitting
#   and plotting pipelines.
#
# WHY THIS FILE EXISTS
#   Multiple scripts need the same preprocessing steps. Centralizing them here
#   avoids subtle inconsistencies (e.g., factor levels or dose mappings).
# -----------------------------------------------------------------------------

library(dplyr) # Data manipulation verbs (mutate, factor handling, etc.).
library(here)  # Project-root-relative file paths.

#' Load and standardize the processed fitting data.
#'
#' This function reads the processed CSV and performs the minimal transformation
#' steps that every analysis script expects:
#'   1) Scenario and Quantity are converted to ordered factors.
#'   2) Dose is mapped from Scenario (0/10/100 mg/kg).
#'   3) xvals is created as a duplicate of Day (used by the fit functions).
#'
#' @param data_path Optional path to the processed data CSV. When NULL, uses the
#'   project default at data/processed-data/processeddata.csv.
#' @return A list with:
#'   - fitdata: standardized data.frame
#'   - scenarios: ordered Scenario levels
#'   - doses: numeric doses in ascending order
#'
load_fit_data <- function(data_path = NULL) {
  # If no path is provided, fall back to the canonical processed CSV.
  if (is.null(data_path)) {
    data_path <- here::here("data", "processed-data", "processeddata.csv")
  }

  # Read the CSV as a plain data.frame; we control factors explicitly below.
  fitdata <- read.csv(data_path, stringsAsFactors = FALSE)

  # Explicit factor level ordering guarantees consistent plotting and fitting.
  scenario_levels <- c("NoTreatment", "PanCytoVir10mg", "PanCytoVir100mg")
  quantity_levels <- c("LogVirusLoad", "IL6", "WeightLossPerc")

  # Standardize columns used by every downstream script.
  fitdata <- fitdata %>%
    mutate(
      Scenario = factor(Scenario, levels = scenario_levels), # Ordered factor.
      Quantity = factor(Quantity, levels = quantity_levels), # Ordered factor.
      Dose = c(0, 10, 100)[as.numeric(Scenario)], # Map scenario -> dose.
      xvals = Day # Alias required by the fit functions.
    )

  # Return both the data and the commonly-used scenario/dose metadata.
  list(
    fitdata = fitdata,
    scenarios = levels(fitdata$Scenario),
    doses = sort(unique(fitdata$Dose))
  )
}
