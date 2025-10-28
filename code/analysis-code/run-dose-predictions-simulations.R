#############################################################################
# script to run dose prediction simulations for all posterior samples
#############################################################################
# load packages needed by this script and the functions it calls
library(here)
library(future.apply) #to do fits in parallel
library(deSolve)
library(caTools)
library(dplyr)

source(here("code/analysis-code/dose-predictions-simulator-function.R"))


bestfit_list <- readRDS(here("results", "output", "bestfit.Rds"))
nsamp <- length(bestfit_list)
simres_list = vector("list", nsamp)

# select doses for which full time-series trajectories will be recorded
timeseries_doses <- c(1, 10, 100, 1000)


# Parallel plan
# don't need a sequential plan here since I don't need diagnostics
# can always initialize more workers and if nsamp is smaller, only a few are used
workers <- 25
plan(multisession, workers = workers)


simres_list <- future_lapply(
  seq_len(nsamp),
  function(i) {
    message(sprintf("processing sample %d", i)) # logs may print slightly out of order
    bestfit <- bestfit_list[[i]]
    simulate_dose_predictions(bestfit, ts_doses = timeseries_doses)
  },
  future.seed = TRUE
)

# switch back to sequential when done
plan(sequential)


# attach the selected doses as metadata for downstream scripts
attr(simres_list, "ts_doses") <- timeseries_doses


saveRDS(
  simres_list,
  here::here('results', 'output', 'dose-response-results.Rds')
)
