# -----------------------------------------------------------------------------
# run-dose-predictions.R
# -----------------------------------------------------------------------------
# PURPOSE
#   Simulate dose-response trajectories using fitted model parameters.
#
# OUTPUTS
#   - results/output/<model>-dose-response-results.Rds
#     This file contains summary AUCs and full time-series for selected doses.
# -----------------------------------------------------------------------------
#
# DETAILED WALKTHROUGH
#   - "User settings": select model, solver settings, and doses to save.
#   - "Load bestfit list": read fitted parameters from disk.
#   - "Simulator selection": pick the model-specific simulator function.
#   - "Run simulations": simulate each bestfit (parallel if multiple).
#   - "Save output": write a single RDS used by plotting scripts.

rm(list = ls(all.names = TRUE)) # Clear workspace for a clean run.

library(here)         # Project-root-relative paths.
library(deSolve)      # ODE solver interface (used by the simulator).
library(caTools)      # trapz() for AUC calculations.
library(dplyr)        # Data manipulation helpers.
library(future)       # Parallel plan configuration.
library(future.apply) # Parallel lapply.

# Model simulators
source(here::here("code", "analysis-code", "functions", "model1-simulator-function.R"))
source(here::here("code", "analysis-code", "functions", "model2-simulator-function.R"))
# Unified fit function (required by build_model_config)
source(here::here("code", "analysis-code", "functions", "fit-function.R"))

# Shared helpers
source(here::here("code", "analysis-code", "functions", "model-config-function.R"))
source(here::here("code", "analysis-code", "functions", "dose-predictions-function.R"))

# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------
model_choice <- "model1" # "model1" or "model2"

# Use the standardized bestfit-sample output (first element is base fit).
bestfit_file <- here::here("results", "output", paste0(model_choice, "-bestfit-sample.Rds"))

# Doses for which full time-series trajectories are retained.
timeseries_doses <- c(0, 1, 10, 1e2, 1e3, 1e4)

# Parallel workers.
workers <- NULL

# ODE solver settings.
solvertype <- "lsoda"
tols <- 1e-9
dt <- 0.01
tfinal <- 7

# -----------------------------------------------------------------------------
# Load bestfit list
# -----------------------------------------------------------------------------
start_time <- Sys.time()
message("Starting run-dose-predictions for model: ", model_choice)
if (!file.exists(bestfit_file)) {
  stop("Bestfit file not found: ", bestfit_file)
}

# Load the list of bestfit objects (one per sample).
bestfit_list <- readRDS(bestfit_file)
nsamp <- length(bestfit_list)

# -----------------------------------------------------------------------------
# Get simulator function from model configuration
# -----------------------------------------------------------------------------
config <- build_model_config(model_choice)
simulatorname <- config$simulatorname

# -----------------------------------------------------------------------------
# Run simulations (parallel if multiple samples)
# -----------------------------------------------------------------------------
if (nsamp > 1) {
  if (is.null(workers)) {
    workers <- max(1, future::availableCores() - 1)
  }
  plan(multisession, workers = workers) # Parallel evaluation.

  simres_list <- future_lapply(
    seq_len(nsamp),
    function(i) {
      message(sprintf("processing sample %d", i))
      bestfit <- bestfit_list[[i]]
      simulate_dose_predictions(
        bestfit,
        ts_doses = timeseries_doses,
        solvertype = solvertype,
        tols = tols,
        dt = dt,
        tfinal = tfinal,
        simulatorname = simulatorname
      )
    },
    future.seed = TRUE
  )

  plan(sequential) # Always reset to sequential.
} else {
  simres_list <- vector("list", 1)
  simres_list[[1]] <- simulate_dose_predictions(
    bestfit_list[[1]],
    ts_doses = timeseries_doses,
    solvertype = solvertype,
    tols = tols,
    dt = dt,
    tfinal = tfinal,
    simulatorname = simulatorname
  )
}

# Attach the selected time-series doses as metadata.
attr(simres_list, "ts_doses") <- timeseries_doses
attr(simres_list, "bestfit_file") <- bestfit_file

# Ensure output directory exists.
output_dir <- here::here("results", "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

output_file <- here::here(
  "results",
  "output",
  paste0(model_choice, "-dose-response-results.Rds")
)

# Save outputs for downstream plotting scripts.
saveRDS(simres_list, output_file)
message("Saved dose-response results to: ", output_file)
end_time <- Sys.time()
elapsed_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
message("run-dose-predictions finished for model: ", model_choice)
message(sprintf("Total elapsed time: %.2f minutes", elapsed_minutes))
