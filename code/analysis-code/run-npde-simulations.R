###############################################################################
# run-npde-simulations.R
#
# This script executes the complete Normalized Prediction Distribution Error
# (NPDE) workflow for the probenecid QSP model.  The goal is to evaluate how
# well the posterior predictive distribution reproduces the experimental data.
#
# The analysis proceeds in four main steps:
#   1. Load the posterior parameter samples (`bestfit.Rds`) and the processed
#      experimental measurements.
#   2. Run the deterministic simulator for every posterior sample and treatment
#      scenario to obtain model predictions at the observation times.
#   3. Combine the simulated predictions with an estimate of measurement noise
#      to compute NPDE values for each individual observation.
#   4. Save the detailed NPDE table, summary diagnostics, and supporting
#      information so that downstream scripts (for example, plotting scripts)
#      can use them without re-running the simulations.
#
# Throughout the script we favour explicit, beginner-friendly code and include
# thorough explanations of the statistical ideas involved in NPDE analysis.
###############################################################################

############################################
## 1. Load required packages and functions
############################################

library(here)    # reliable file paths that work across operating systems
library(dplyr)   # data wrangling verbs (mutate, summarise, join, ...)
library(tidyr)   # reshaping data (pivot_longer, etc.)
library(purrr)   # functional programming helpers (map, imap)
library(readr)   # reading CSV files
library(tibble)  # building tidy tibbles and list-columns

# Core simulator and the NPDE helper functions defined for this project
source(here::here("code/analysis-code/model-simulator-function.R"))
source(here::here("code/analysis-code/npde-simulation-function.R"))

############################################
## 2. Read posterior samples and data
############################################

bestfit_path <- here::here("results", "output", "bestfit.Rds")
if (!file.exists(bestfit_path)) {
  stop("The file 'results/output/bestfit.Rds' was not found. Run the model fitting workflow before computing NPDE diagnostics.")
}

bestfit_list <- readRDS(bestfit_path)
if (length(bestfit_list) == 0) {
  stop("No posterior samples found in 'bestfit.Rds'. NPDE calculations require at least one sample.")
}

# Load the processed experimental data using the helper function.  The function
# provides additional convenience columns (scenario ordering, dose values, ...)
# that the rest of the script relies on.
data_path <- here::here("data", "processed-data", "processeddata.csv")
observed_data <- load_observed_dataset(data_path)

############################################
## 3. Run the NPDE analysis pipeline
############################################

npde_results <- run_full_npde_analysis(
  bestfit_list = bestfit_list,
  data_path = data_path,
  dt = 0.01
)

############################################
## 4. Persist results to disk for later use
############################################

output_dir <- here::here("results", "output")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

tables_dir <- here::here("results", "tables")
if (!dir.exists(tables_dir)) {
  dir.create(tables_dir, recursive = TRUE)
}

# Save the full results as an RDS file.  This preserves the nested list
# structure, including intermediate objects such as the posterior predictions.
saveRDS(
  npde_results,
  file = file.path(output_dir, "npde-results.Rds")
)

# Export the detailed NPDE table and the summary diagnostics as CSV files so
# they can be inspected in spreadsheets or external tools.
readr::write_csv(
  npde_results$npde_table,
  file = file.path(tables_dir, "npde-values.csv")
)

readr::write_csv(
  npde_results$summary,
  file = file.path(tables_dir, "npde-summary.csv")
)

############################################
## 5. User-friendly console output
############################################

message("NPDE analysis complete.")
message("- Saved detailed results to results/output/npde-results.Rds")
message("- Saved per-observation NPDE values to results/tables/npde-values.csv")
message("- Saved summary diagnostics to results/tables/npde-summary.csv")

###############################################################################
# End of script
###############################################################################
