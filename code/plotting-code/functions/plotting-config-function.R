# -----------------------------------------------------------------------------
# plotting-config-function.R
# -----------------------------------------------------------------------------
# PURPOSE
#   Provide shared labels, colors, and factor mappings for plotting. Keeping
#   these values in one place avoids duplicated hard-coded strings and reduces
#   the risk of inconsistent labels across figures.
# -----------------------------------------------------------------------------

#' Return shared plotting configuration values.
#'
#' @return A list containing:
#'   - scenario_map: mapping of numeric dose to scenario labels
#'   - scenario_levels: ordered scenario labels
#'   - quantity_levels: ordered measurement labels
#'   - var_colors: colors for variables
#'   - var_labs: display labels for variables
#'   - scen_labs: display labels for scenarios
#'
get_plotting_config <- function() {
  # Map numeric dose values to scenario labels used in the dataset.
  scenario_map <- c(
    `0` = "NoTreatment",
    `10` = "PanCytoVir10mg",
    `100` = "PanCytoVir100mg"
  )

  # Return a named list so downstream scripts can pick what they need.
  list(
    scenario_map = scenario_map,
    scenario_levels = c("NoTreatment", "PanCytoVir10mg", "PanCytoVir100mg"),
    quantity_levels = c("LogVirusLoad", "IL6", "WeightLossPerc"),
    var_colors = c(
      LogVirusLoad = "#0072B2",
      IL6 = "#D55E00",
      WeightLossPerc = "#009E73"
    ),
    var_labs = c(
      LogVirusLoad = "Log Virus Load",
      IL6 = "IL-6",
      WeightLossPerc = "Weight Loss (%)"
    ),
    scen_labs = c(
      NoTreatment = "No Treatment",
      PanCytoVir10mg = "10 mg/kg",
      PanCytoVir100mg = "100 mg/kg"
    )
  )
}
