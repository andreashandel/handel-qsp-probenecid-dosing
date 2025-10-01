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


# Parallel plan (Windows-safe)
#workers <- max(1, parallel::detectCores(logical = TRUE) - 1)
workers <- 25
plan(multisession, workers = workers)


simres_list <- future_lapply(
  seq_len(nsamp),
  function(i) {
    message(sprintf("processing sample %d", i)) # logs may print slightly out of order
    bestfit <- bestfit_list[[i]]
    simulate_dose_predictions(bestfit)
  },
  future.seed = TRUE
)

# switch back to sequential when done
plan(sequential)


saveRDS(
  simres_list,
  here::here('results', 'output', 'dose-response-results.Rds')
)
