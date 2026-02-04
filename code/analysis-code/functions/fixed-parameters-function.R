# -----------------------------------------------------------------------------
# fixed-parameters-function.R
# -----------------------------------------------------------------------------
# PURPOSE
#   Load the fixed (non-fitted) parameter values and their descriptive labels
#   from CSV files stored under data/processed-data/.
#
# WHY THIS FILE EXISTS
#   Multiple scripts need the same CSV parsing and name handling. Centralizing
#   this logic prevents inconsistent naming or factor handling across scripts.
# -----------------------------------------------------------------------------

library(here) # Project-root-relative file paths.

#' Load fixed parameters from a CSV file.
#'
#' @param csv_path Path to the fixed-parameter CSV file. If NULL, this function
#'   will error because it does not know which model to load.
#' @return A list with:
#'   - values: named numeric vector of fixed parameter values
#'   - labels: named character vector of descriptive labels
#'   - raw: the original data.frame (useful for debugging/inspection)
#'
load_fixed_parameters <- function(csv_path) {
  # Guardrail: this function needs a file path to operate.
  if (is.null(csv_path)) {
    stop("csv_path must be provided.")
  }

  # Read the CSV without implicit factors for full control of types.
  fixedpars_df <- read.csv(csv_path, stringsAsFactors = FALSE)

  # Normalize whitespace in parameter names and labels.
  fixedpars_df$parname <- trimws(fixedpars_df$parname)
  fixedpars_df$parnamefull <- trimws(fixedpars_df$parnamefull)

  # Build a named numeric vector of parameter values.
  fixed_values <- as.numeric(fixedpars_df$value)
  names(fixed_values) <- fixedpars_df$parname

  # Build a named vector of human-readable labels for UI/table output.
  fixed_labels <- fixedpars_df$parnamefull
  names(fixed_labels) <- fixedpars_df$parname

  # Return values, labels, and the original data.frame for inspection.
  list(
    values = fixed_values,
    labels = fixed_labels,
    raw = fixedpars_df
  )
}
