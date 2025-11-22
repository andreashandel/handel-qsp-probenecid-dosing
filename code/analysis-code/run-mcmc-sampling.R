#############################################################################
# Bayesian parameter sampling via BayesianTools
#############################################################################
# Goal: leverage the existing likelihood/ODE model to obtain posterior samples
# for fitted parameters using an off-the-shelf MCMC implementation. We wrap the
# current log-likelihood in a BayesianTools "Bayesian setup" and run the DREAMzs
# sampler, which adapts automatically and runs multiple chains in parallel.
# This avoids maintaining a custom sampler while fully reusing existing code.
#############################################################################

# ---------------------------------------------------------------------------
# Required packages
# ---------------------------------------------------------------------------
# here          : consistent file paths
# dplyr         : light data wrangling while loading processed data
# deSolve       : needed because likelihood evaluates ODEs
# BayesianTools : provides generic MCMC samplers (DE-MCMC, DREAMzs, etc.)
# future        : optional; only used to parallelize chains via BayesianTools
# ---------------------------------------------------------------------------
library(here)
library(dplyr)
library(deSolve)
library(BayesianTools)

# ---------------------------------------------------------------------------
# Source underlying simulator and likelihood code so the Bayesian workflow
# matches the maximum-likelihood implementation exactly.
# ---------------------------------------------------------------------------
source(here("code", "analysis-code", "model-simulator-function.R"))
source(here("code", "analysis-code", "fit-model-function.R"))

# Load processed data and parameter meta-data --------------------------------
fitdata <- read.csv(here("data", "processed-data", "processeddata.csv"))
fitdata$Scenario <- factor(
  fitdata$Scenario,
  levels = c("NoTreatment", "PanCytoVir10mg", "PanCytoVir100mg")
)
fitdata$Quantity <- factor(
  fitdata$Quantity,
  levels = c("LogVirusLoad", "IL6", "WeightLossPerc")
)
fitdata$Dose <- c(0, 10, 100)[as.numeric(fitdata$Scenario)]
fitdata$xvals <- fitdata$Day

fixedparsdata <- read.csv(here("data", "processed-data", "fixed-parameters.csv"))
bestfit_list <- readRDS(here("results", "output", "bestfit.Rds"))
bestfit_baseline <- bestfit_list[[1]]

Y0 <- bestfit_baseline$Y0
fixedpars <- bestfit_baseline$fixedpars
fitparnames <- bestfit_baseline$fitparnames
par_ini <- as.numeric(bestfit_baseline$solution)
names(par_ini) <- fitparnames
doses <- unique(fitdata$Dose)
scenarios <- levels(fitdata$Scenario)
solvertype <- "vode"
tols <- 1e-9
tfinal <- 7
dt <- 0.02

# ---------------------------------------------------------------------------
# Parameter bounds (uniform priors) – identical to run-model-fitting.R
# ---------------------------------------------------------------------------
b = 1e-8
bl = 1e-12
bh = 1e-2
k = 1e-4
kl = 1e-8
kh = 1e5
p = 2e3
pl = 1e-3
ph = 1e10
kF = 1
kFl = 1e-10
kFh = 1e3
cV = 10
cVl = 0.01
cVh = 1e5
gF = 0.1
gFl = 1e-3
gFh = 1e3
hV = 1e4
hVl = 1e-5
hVh = 1e8
Fmax = 5
Fmaxl = 0.1
Fmaxh = 1000
hF = 1
hFl = 1e-3
hFh = 1e5
gS = 10
gSl = 1e-4
gSh = 1e4
cS = 1
cSl = 1e-2
cSh = 1e4
Emax_F = 0.5
Emax_Fl = 1e-3
Emax_Fh = 1
C50_F = 1
C50_Fl = 1e-7
C50_Fh = 1e6
C50_V = 1
C50_Vl = 1e-7
C50_Vh = 1e6

var_by_qty <- fitdata %>%
  group_by(Quantity) %>%
  summarize(v = var(Value, na.rm = TRUE), .groups = "drop")
sigma_all <- c(
  sigma_add_LogVirusLoad = sqrt(as.numeric(var_by_qty[1, 2])),
  sigma_prop_LogVirusLoad = 0.0,
  sigma_add_IL6 = sqrt(as.numeric(var_by_qty[2, 2])),
  sigma_prop_IL6 = 0.0,
  sigma_add_WeightLossPerc = sqrt(as.numeric(var_by_qty[3, 2])),
  sigma_prop_WeightLossPerc = 0.0
)
sigma_to_fit <- c(
  "sigma_add_LogVirusLoad",
  "sigma_add_IL6",
  "sigma_add_WeightLossPerc"
)
sigma_fit_ini <- sigma_all[sigma_to_fit]

lb <- c(
  bl, kl, pl, kFl, cVl, gFl, hVl, Fmaxl, hFl, gSl, cSl, Emax_Fl, C50_Fl, C50_Vl,
  rep(1e-6, length(sigma_fit_ini))
)
ub <- c(
  bh, kh, ph, kFh, cVh, gFh, hVh, Fmaxh, hFh, gSh, cSh, Emax_Fh, C50_Fh, C50_Vh,
  rep(1e3, length(sigma_fit_ini))
)
names(lb) <- names(ub) <- fitparnames

# ---------------------------------------------------------------------------
# Log-likelihood wrapper
# ---------------------------------------------------------------------------
# BayesianTools expects a function mapping a parameter vector to log-likelihood.
# We simply call fit_model_function() and return its negative (since the fitter
# returns an objective to minimize).
# ---------------------------------------------------------------------------
loglik_fun <- function(theta) {
  theta <- as.numeric(theta)
  names(theta) <- fitparnames
  if (any(theta <= lb) || any(theta >= ub)) {
    return(-Inf)
  }
  obj <- fit_model_function(
    params = theta,
    fitdata = fitdata,
    Y0 = Y0,
    tfinal = tfinal,
    dt = dt,
    fitparnames = fitparnames,
    fixedpars = fixedpars,
    doses = doses,
    scenarios = scenarios,
    solvertype = solvertype,
    tols = tols
  )
  return(-obj)
}

# Uniform prior across parameter bounds -------------------------------------
prior_fun <- createUniformPrior(lower = lb, upper = ub, best = par_ini)

# Combine into BayesianTools setup ------------------------------------------
bt_setup <- createBayesianSetup(
  likelihood = loglik_fun,
  prior = prior_fun,
  names = fitparnames
)

# MCMC settings -------------------------------------------------------------
mcmc_settings <- list(
  iterations = 20000,
  nrChains = 4,
  message = TRUE,
  thin = 1
)
set.seed(20250306)

# Run DREAMzs sampler (good default for correlated posteriors) ---------------
mcmc_out <- runMCMC(
  bayesianSetup = bt_setup,
  sampler = "DREAMzs",
  settings = mcmc_settings
)

# Posterior summary ---------------------------------------------------------
summary_out <- summary(mcmc_out)
print(summary_out)

# Save chains and summary for downstream use -------------------------------
saveRDS(
  list(
    sampler_output = mcmc_out,
    summary = summary_out,
    settings = mcmc_settings
  ),
  here("results", "output", "mcmc-samples.Rds")
)

cat("BayesianTools sampling complete.\n")
