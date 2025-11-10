# R script to fit data to the QSP model

# clear workspace
rm(list = ls(all.names = TRUE))


# various packages needed by this script and the function it calls
library(here)
library(dplyr)
library(nloptr)
library(deSolve)
library(lhs) # for Latin Hypercube Sampling of fixed parameters
library(future)
library(future.apply) #to do fits in parallel
library(beepr) # to make a sound when fitting is done

# this file contains the ODE model as a function
source(here::here('code/analysis-code/model-simulator-function.R'))

# function that performs a single fit iteration
source(here::here('code/analysis-code/fit-model-function.R'))

#load and process data
file_path = here::here("data/processed-data/processeddata.csv")
fitdata = read.csv(file_path)

# load fixed parameters from file
# this does not include potentially fixed parameters for error distributions (i.e. the sigmas)
pars_file = here::here('data/processed-data/fixed-parameters.csv')
fixedparsdata = read.csv(pars_file)


# make Scenario an ordered factor
fitdata$Scenario = factor(
  fitdata$Scenario,
  levels = c("NoTreatment", "PanCytoVir10mg", "PanCytoVir100mg")
)
# the 3 different treatment scenarios
scenarios = unique(fitdata$Scenario)

# make sure the different variables are alwaysi in the same oder
fitdata$Quantity = factor(
  fitdata$Quantity,
  levels = c("LogVirusLoad", "IL6", "WeightLossPerc")
)


# add Dose to fitdata dataframe
# values should be 0, 10 and 100 for the three scenarios
# this is in mg/kg
# scaling to actual amount happens inside simulator 
fitdata$Dose = c(0, 10, 100)[as.numeric(fitdata$Scenario)]

# rename Day to xvals (by adding a column called xvals)
fitdata$xvals = fitdata$Day


# starting values for variables
U = 1e7 #from 2018 PCB paper
I = 0
V = 1 #assuming 1 virion to start
F = 0 #no innate initially
A = 0 #initial number of activated adaptive response
S = 0 #no symptoms
Ad = 0 # no drug
Ac = 0 #no drug
At = 0 #no drug

#combine initial conditions into a vector
Y0 = c(Ad = Ad, Ac = Ac, At = At, U = U, I = I, V = V, F = F, A = A, S = S)

#defining parameters to fit
########################################
# starting values for fitted parameters
# also defining low and high/upper bounds
# see manuscript for their definitions
########################################
b = 1e-8 # infection rate
bl = 1e-12
bh = 1e-2
k = 1e-4 #adaptive response virus removal
kl = 1e-8
kh = 1e5
p = 2e3 #virus production rate
pl = 1e-3
ph = 1e10
kF = 1 #innate impact on virus production
kFl = 1e-10
kFh = 1e3
cV = 10 #virus clearance rate
cVl = 0.01
cVh = 1e5

gF = 0.1 #max innate growth
gFl = 1e-3
gFh = 1e3
hV = 1e4 # saturation for virus induction effect
# hVl = 1e-5
# hVh = 1e8
Fmax = 5 # max innate response
Fmaxl = 0.1
Fmaxh = 10
hF = 10 # T-cell induction response
hFl = 100
hFh = 1e5
gS = 10 #induction of symptoms by innate
gSl = 1e-4
gSh = 1e4
cS = 1 # rate of symptom decline
cSl = 1e-2
cSh = 1e4

#PD
Emax_F = 0.5 #strength of reduction of innate response by drug
Emax_Fl = 1e-3
Emax_Fh = 1
C50_F = 1 #50% reduction on innnate
C50_Fl = 1e-7
C50_Fh = 1e6
C50_V = 1 #50% reduction on virus
C50_Vl = 1e-7
C50_Vh = 1e6

# combine all main fitted parameters into a vector
par_ini_full = c(
  b = b,
  k = k,
  p = p,
  kF = kF,
  cV = cV,
  gF = gF,
  Fmax = Fmax,
  hF = hF,
  gS = gS,
  cS = cS,
  Emax_F = Emax_F,
  C50_F = C50_F,
  C50_V = C50_V
)

# Compute empirical variance per Quantity; small floor to avoid zero
var_by_qty <- fitdata %>%
  group_by(Quantity) %>%
  summarize(v = var(Value, na.rm = TRUE), .groups = "drop")


### NEW/CHANGED: define ALL six sigma params here (initial values)
sigma_all <- c(
  sigma_add_LogVirusLoad = sqrt(as.numeric(var_by_qty[1, 2])),
  sigma_prop_LogVirusLoad = 0.0,
  sigma_add_IL6 = sqrt(as.numeric(var_by_qty[2, 2])),
  sigma_prop_IL6 = 0.0,
  sigma_add_WeightLossPerc = as.numeric(sqrt(var_by_qty[3, 2])),
  sigma_prop_WeightLossPerc = 0.0
)

### NEW/CHANGED: choose which of the six to FIT (edit this vector as needed)
sigma_to_fit <- c(
  "sigma_add_LogVirusLoad",
  "sigma_add_IL6",
  "sigma_add_WeightLossPerc"
  # e.g., add "sigma_prop_LogVirusLoad" here if you want to fit it, too
)


### NEW/CHANGED: split fitted vs fixed sigmas
sigma_fit_ini <- sigma_all[sigma_to_fit]
sigma_fixed <- sigma_all[setdiff(names(sigma_all), sigma_to_fit)]

### add the fitted sigmas to the fitted parameter vector
par_ini_full <- c(par_ini_full, sigma_fit_ini)

par_ini <- as.numeric(par_ini_full)
fitparnames <- names(par_ini_full)


# this is saved for later production of table
parlabels = c(
  "Virus infection rate",
  "Adaptive response clearance rate",
  "Virus production rate",
  "Innate response supression strength",
  "Virus removal rate",
  "Maximum innate response induction",
  "Maximum innate response",
  "Adaptive response half-maximum induction",
  "Symptom induction rate",
  "Symptom decay rate",
  "Maximum innate response supression",
  "Half maximum of innate response effect",
  "Half maximum of virus suppression effect",
  "Sigma of LogVirusLoad",
  "Simga of IL6",
  "Sigma of WeightLossPerc"
)


# starting values either from best fit values of previous run or values above
# can be commented out if one wants to start
# with the above values
oldbestfit = readRDS(here::here('results', 'output', 'bestfit.Rds'))
par_ini = as.numeric(oldbestfit[[1]]$solution)
# names(par_ini) = oldbestfit[[1]]$fitparnames
#par_ini = c(6.38564317479985e-9, 0.0008139662918231, 45632.927287892, 0.000524055968325585, 126.02633258623, 0.0612491860382419, 5, 100, 69.3937856059817, 0.0100000000000009, 0.00172438718255394, 1.00135950945941e-07, 3.9452670546606e-07, 4.78095447065562, 0.307283559709252, 8.78609676438555) 
  

# upper and lower bounds of parameters
lb = as.numeric(c(
  bl,
  kl,
  pl,
  kFl,
  cVl,
  gFl,
  Fmaxl,
  hFl,
  gSl,
  cSl,
  Emax_Fl,
  C50_Fl,
  C50_Vl
))
ub = as.numeric(c(
  bh,
  kh,
  ph,
  kFh,
  cVh,
  gFh,
  Fmaxh,
  hFh,
  gSh,
  cSh,
  Emax_Fh,
  C50_Fh,
  C50_Vh
))

### append bounds for the sigmas
lb <- c(lb, rep(1e-6, length(sigma_fit_ini)))
ub <- c(ub, rep(1e3, length(sigma_fit_ini)))


# check bounds. if initial conditions are outside bounds, give a warning and adjust to bound.
if (sum((ub - par_ini) < 0) > 0) {
  print(
    "Warning: initial value is larger than upper bound, setting it to upper bound"
  )
  #par_ini = pmin(ub, par_ini)
}
if (sum((par_ini - lb) < 0) > 0) {
  print(
    "Warning: initial value is smaller than lower bound, setting it to lower bound"
  )
  #par_ini = pmax(lb, par_ini)
}

# number of samples
nsamp = 0 # if this is 0, we only fit for the baseline values of the fixed parameters

# settings for optimizer
algname = "NLOPT_LN_COBYLA"
#algname = "NLOPT_LN_NELDERMEAD"
#algname = "NLOPT_LN_SBPLX"
maxsteps = 500 #number of steps/iterations for algorithm
maxtime = 10 * 60 * 60 #maximum time in seconds (h*m*s)
ftol_rel = 1e-10

# settings for ODE solver
solvertype = "vode"
tols = 1e-9
tfinal = 7 #time of last data point
dt = 0.02 # time step for which we want results returned

##############
# make samples for fixed parameters and re-fit for each sample
# this excludes the sigmas
##############

# extract relevant content from fixed parameters object
# the object was loaded at the beginning of the script
fixedpars = fixedparsdata[, 3]
names(fixedpars) = fixedparsdata[, 1]


# always add baseline as sample to fit
samples_list <- list(fixedpars)

# create a list of samples for fixed parameters ranging from half to double, uniformly sampled with LHS
set.seed(1234) #for reproducibility of sampling and also for parallel fitting (which should not require randomness)
if (nsamp > 0) {
  lower <- 0.5 * fixedpars
  upper <- 2 * fixedpars
  U <- lhs::randomLHS(nsamp, length(fixedpars)) # nsamp x p in [0,1]
  M <- sweep(U, 2, (upper - lower), "*")
  M <- sweep(M, 2, lower, "+")
  colnames(M) <- names(fixedpars)
  samples_list_2 <- lapply(seq_len(nsamp), function(i) {
    setNames(as.numeric(M[i, ]), names(fixedpars))
  })
  samples_list = c(samples_list, samples_list_2)
}

# Emax_V is not sampled, set back to 1 everywhere
samples_list <- lapply(samples_list, function(x) {
  x["Emax_V"] <- 1
  x
})

### NEW/CHANGED: append the FIXED sigmas (NOT sampled) to each fixed-pars sample
samples_list <- lapply(samples_list, function(x) c(x, sigma_fixed))

# make an empty list of length nsamp + 1 (always fit baseline)
bestfit_all = vector("list", nsamp + 1)

# record current time
start_time <- proc.time()


#################################################
# --- function to process a single sample i ---
# needed for parallel processing
#################################################
eval_one_sample <- function(i, print_level) {
  message(sprintf(
    "processing sample %d at %s",
    i,
    format(Sys.time(), "%H:%M:%S")
  ))

  fixedpars_i <- samples_list[[i]]

  # ---- fit ----
  bestfit <- nloptr::nloptr( x0 = par_ini,
    eval_f = fit_model_function, lb = lb, ub = ub,
    opts = list(
      algorithm = algname,
      maxeval = maxsteps,
      print_level = print_level,
      ftol_rel = ftol_rel
    ),
    fitdata = fitdata,
    Y0 = Y0,
    tfinal = tfinal,
    dt = dt,
    fitparnames = fitparnames,
    fixedpars = fixedpars_i,
    doses = unique(fitdata$Dose),
    scenarios = scenarios,
    solvertype = solvertype,
    tols = tols
  )
  # finished with fitting
  # doing some after fitting stuff

  # extract params
  params <- bestfit$solution
  names(params) <- fitparnames

  # pack extras
  parstring <- paste0("c(", paste(as.numeric(params), collapse = ", "), ")")
  bestfit$parstring <- parstring #same as solution but formatted so we can stick i in easily as start value
  bestfit$fitpars <- params #same as solution but with names for parameters
  bestfit$fitparnames <- fitparnames
  bestfit$fixedpars <- fixedpars_i
  bestfit$Y0 <- Y0
  bestfit$fitdata <- fitdata
  bestfit$parlabels <- parlabels #full names/labels for parameters
  bestfit$algorithm <- algname

  return(bestfit) # return the finished best fit
}


#########################################
# do things either in parallel or not
#########################################
if (length(samples_list) > 1) {
  n_workers <- 25
  workers <- min(n_workers, future::availableCores())
  print_level <- 0
  future::plan(multisession, workers = workers)
  message("Running in parallel with ", workers, " workers.")

  # --- run in parallel; reproducible RNG across workers ---
  bestfit_all <- future_lapply(
    seq_along(samples_list),
    eval_one_sample,
    print_level,
    future.seed = TRUE
  )

  # optional: switch back to sequential when done
  future::plan(sequential)
} else {
  message("Single sample detected; running sequentially (no futures).")
  print_level <- 1 #diagnostics

  bestfit_all[[1]] <- eval_one_sample(1, print_level)
}


# record end time and print total time elapsed in minutes
#capture time taken for fit
tdiff = proc.time() - start_time
runtime_minutes = tdiff[[3]] / 60
cat('model fit took this many minutes:', runtime_minutes, '\n')
cat('************** \n')
cat('used algorithm: ', algname, '\n')
cat('************** \n')
# Safe print (won’t error if oldbestfit wasn’t loaded)
if (exists("oldbestfit")) {
  cat('initial objective function: ', oldbestfit[[1]]$objective, '\n')
} else {
  cat('initial objective function: (no previous run loaded)\n')
}
cat('************** \n')
cat('final objective functions: ', '\n')
# print all objective function values for each sample
for (i in 1:(nsamp + 1)) {
  print(bestfit_all[[i]]$objective)
}
# play a sound when done
beepr::beep(2)
# statement that prints best fit parameter values as vector
print(bestfit_all[[1]]$parstring)

# copy prior best fit to a new file if it exists
if (file.exists(here::here('results', 'output', 'bestfit.Rds'))) {
  file.copy(
    from = here::here('results', 'output', 'bestfit.Rds'),
    to = here::here('results', 'output', 'oldbestfit.Rds'),
    overwrite = TRUE
  )
}
# save new best fit
saveRDS(bestfit_all, here::here('results', 'output', 'bestfit.Rds'))
