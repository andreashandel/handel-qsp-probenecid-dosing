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
model_choice <- "model2" # "model1" or "model2"

# Use the standardized bestfit-sample output (first element is base fit).
bestfit_file <- here::here("results", "output", paste0(model_choice, "-bestfit-sample.Rds"))

#x <- readRDS(here::here("results", "output", paste0(model_choice, "-bestfit-multistart.Rds")))

#Remove fixed parameters from the app UI. And make sure that the app treats parameters on the right scale. There is no reason to ever do a log transform of the parameters for the app, so keep all on a linear scale and make sure the functions which are called also treat the model parameters on the linear scale. The objective value is still not aligned. The value reported in the UI right now is 9.38717e+01 while the 

# Doses for which full time-series trajectories are retained.
timeseries_doses <- c(0, 1, 10, 1e2, 1e3, 1e4)

# All doses to simulate (must include timeseries_doses).
dose_grid <- 10^seq(-2, 4, length = 50)
all_doses <- sort(unique(c(timeseries_doses, dose_grid)))

# Dosing schedules.
schedule_defs <- list(
  s1 = list(txstart = 1, txend = 3.9, txinterval = 0.5, name = "s1", label = "baseline"),
  s2 = list(txstart = 2, txend = 4.9, txinterval = 0.5, name = "s2", label = "d2 start"),
  s3 = list(txstart = 3, txend = 5.9, txinterval = 0.5, name = "s3", label = "d3 start"),
  s4 = list(txstart = 1, txend = 3.9, txinterval = 1, name = "s4", label = "daily tx"),
  s5 = list(txstart = 1, txend = 1, txinterval = 1, name = "s5", label = "single tx")
)

# Parallel workers.
workers <- 25

# ODE solver settings.
solvertype <- "vode"
tols <- 1e-10
dt <- 0.02
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
        all_doses = all_doses,
        schedule_defs = schedule_defs,
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
    all_doses = all_doses,
    schedule_defs = schedule_defs,
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

# -----------------------------------------------------------------------------
# Save warnings/errors from dose-response simulations
# -----------------------------------------------------------------------------
message_rows <- lapply(seq_along(simres_list), function(i) {
  res <- simres_list[[i]]
  warn_df <- res$warnings
  fail_df <- res$failures

  warn_df$Type <- if (nrow(warn_df)) "warning" else character(0)
  fail_df$Type <- if (nrow(fail_df)) "error" else character(0)

  bind_rows(warn_df, fail_df) %>%
    mutate(Sample = i)
})

messages_df <- bind_rows(message_rows)

messages_file <- here::here(
  "results",
  "output",
  paste0(model_choice, "-dose-response-messages.csv")
)

if (nrow(messages_df)) {
  write.csv(messages_df, messages_file, row.names = FALSE)
  message("Saved dose-response warnings/errors to: ", messages_file)
}

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
