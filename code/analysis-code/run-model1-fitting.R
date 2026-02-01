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
source(here::here('code/analysis-code/model1-simulator-function.R'))

# function that performs a single fit iteration
source(here::here('code/analysis-code/model1-fit-function.R'))

#load and process data
file_path = here::here("data/processed-data/processeddata.csv")
fitdata = read.csv(file_path)

# load fixed parameters from file
# this does not include potentially fixed parameters for error distributions (i.e. the sigmas)
pars_file = here::here('data/processed-data/model1-fixed-parameters.csv')
fixedparsdata = read.csv(pars_file)

# number of samples
nsamp = 0 # if this is 0, we only fit for the baseline values of the fixed parameters
n_workers <- 34 #number of workers for parallel processing - is ignored for nsamp = 0


# load prior best fit, can be used as starting condition
if (nsamp == 0) {
  bestfitfile = 'model1-bestfit-single.Rds'
} else {
  bestfitfile = 'model1-bestfit-sample.Rds'
}
# check if prior best fit file exists, then load
if (file.exists(here::here('results', 'output', bestfitfile))) {
    oldbestfit = readRDS(here::here('results', 'output', bestfitfile))
} else {
  oldbestfit = NULL
}


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
bh = 1e-5
k = 1e-5 #adaptive response virus removal
kl = 1e-10
kh = 1
p = 1e4 #virus production rate
pl = 1
ph = 1e7
kF = 0.1 #innate impact on virus production
kFl = 1e-1
kFh = 1e2
cV = 100 #virus clearance rate
cVl = 0.1
cVh = 1e5

gF = 1 #max innate growth
gFl = 1e-3
gFh = 1e3
hV = 1e3 # saturation for virus induction effect
hVl = 1e-2
hVh = 1e5
Fmax = 2 # max innate response
Fmaxl = 0.1
Fmaxh = 1e3
hF = 1 # T-cell induction response
hFl = 1e-5
hFh = 1e3
gS = 10 #induction of symptoms by innate
gSl = 1e-3
gSh = 1e3
cS = 1 # rate of symptom decline
cSl = 1e-3
cSh = 1e3

#PD
Emax_F = 1 #strength of reduction of innate response by drug
Emax_Fl = 1e-3
Emax_Fh = 1
C50_F = 1e-5 #50% reduction on innnate
C50_Fl = 1e-10
C50_Fh = 1e2
C50_V = 1e-8 
C50_Vl = 1e-10
C50_Vh = 1e2


##########################
# Code bits for additive or multiplicative errors
##########################

# Compute empirical variance per Quantity; small floor to avoid zero
var_by_qty <- fitdata %>%
  group_by(Quantity) %>%
  summarize(v = var(Value, na.rm = TRUE), .groups = "drop")


### set initial values for all six sigmas for additive and multiplicative errors params here 
### any zero means that error term is not used

sigma_all <- c(
  sigma_add_LogVirusLoad = sqrt(as.numeric(var_by_qty[1, 2])),
  sigma_prop_LogVirusLoad = 0.0,
  sigma_add_IL6 = sqrt(as.numeric(var_by_qty[2, 2])),
  sigma_prop_IL6 = 0.0,
  sigma_add_WeightLossPerc = as.numeric(sqrt(var_by_qty[3, 2])),
  sigma_prop_WeightLossPerc = 0.0
)

### choose which of the six to FIT (edit this vector as needed)
sigma_to_fit <- c(
  #"sigma_add_LogVirusLoad",
  #"sigma_add_IL6",
  #"sigma_add_WeightLossPerc"
  # e.g., add "sigma_prop_LogVirusLoad" here if you want to fit it, too
)

### split fitted vs fixed sigmas
sigma_fit_ini <- sigma_all[sigma_to_fit]
sigma_fixed <- sigma_all[setdiff(names(sigma_all), sigma_to_fit)]



# combine all main fitted parameters into a vector
par_ini_full = c(
  b = b,
  k = k,
  p = p,
  kF = kF,
  cV = cV,
  gF = gF,
  hV = hV,
  Fmax = Fmax,
  hF = hF,
  gS = gS,
  cS = cS,
  Emax_F = Emax_F,
  C50_F = C50_F,
  C50_V = C50_V
)

# upper and lower bounds of parameters
lb = as.numeric(c(
  bl,
  kl,
  pl,
  kFl,
  cVl,
  gFl,
  hVl,
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
  hVh,
  Fmaxh,
  hFh,
  gSh,
  cSh,
  Emax_Fh,
  C50_Fh,
  C50_Vh
))

### add the fitted sigmas to the fitted parameter vector
par_ini_full <- c(par_ini_full, sigma_fit_ini)
### append bounds for the sigmas
lb <- c(lb, rep(1e-6, length(sigma_fit_ini)))
ub <- c(ub, rep(1e3, length(sigma_fit_ini)))

# give bounds the same names as the fit parameter vector
names(lb) <- names(par_ini_full)
names(ub) <- names(par_ini_full)

# this is saved for later production of table
parlabels = c(
  b = "Virus infection rate",
  k = "Adaptive response clearance rate",
  p = "Virus production rate",
  kF = "Innate response supression strength",
  cV = "Virus removal rate",
  gF = "Maximum innate response induction",
  hV = "Adaptive response half-maximum induction",
  Fmax = "Maximum innate response",
  hF = "Adaptive response half-maximum induction",
  gS = "Symptom induction rate",
  cS = "Symptom decay rate",
  Emax_F = "Maximum drug effect on innate response",
  C50_F = "Half maximum of innate response effect",
  C50_V = "Half maximum of virus suppression effect",
  sigma_add_LogVirusLoad = "Sigma of LogVirusLoad",
  sigma_add_IL6 = "Simga of IL6",
  sigma_add_WeightLossPerc = "Sigma of WeightLossPerc"
)



###############################################################
# optionally hold selected parameters fixed for testing.
# Provide parameter names in `user_fixed_param_names` to exclude
# them from the fitted vector and treat them as additional fixed
# parameters later in the script. Example below fixes hV and C50_F.
###############################################################
#user_fixed_params <- c()
user_fixed_params <- c(Emax_F = 1)
#user_fixed_params <- c(Emax_F = 1, hF = 1, p = 210)
if (length(user_fixed_params)) {
  missing_names <- setdiff(names(user_fixed_params), names(par_ini_full))
  
  par_ini_full <- par_ini_full[setdiff(names(par_ini_full), names(user_fixed_params))]
  # also adjust bounds vectors
  lb <- lb[setdiff(names(par_ini_full), names(user_fixed_params))]
  ub <- ub[setdiff(names(par_ini_full), names(user_fixed_params))]
}

fitparnames <- names(par_ini_full)


# keep labels only for parameters still being fitted
parlabels <- parlabels[fitparnames]
if (any(is.na(parlabels))) {
  stop("Parlabel definitions missing entries for some fitted parameters.")
}

if (length(parlabels) != length(par_ini_full)) {
  stop("length of parlabels does not match length of par_ini")
}



# name of underlying model simulator function
simulatorname = "model1_simulator"



# settings for optimizer
#algname = "NLOPT_LN_COBYLA"
#algname = "NLOPT_LN_BOBYQA"
algname = "NLOPT_LN_NEWUOA"
#algname = "NLOPT_LN_PRAXIS"
#algname = "NLOPT_LN_NELDERMEAD"
#algname = "NLOPT_LN_SBPLX"
maxsteps = 500 #number of steps/iterations for algorithm
maxtime = 10 * 60 * 60 #maximum time in seconds (h*m*s)
ftol_rel = 1e-10
logfit = 1 #if 1 we fit parameters in log space, if 0 in linear space

# settings for ODE solver
#solvertype = "lsoda"
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

# append any user-specified fixed parameters 
if (length(user_fixed_params)) {
  fixedpars <- c(fixedpars, user_fixed_params)
}


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

### append the FIXED sigmas (NOT sampled) to each fixed-pars sample
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

  # starting values either from best fit values of previous run or values above
  # can be commented out if one wants to start
  # with the above values

  # assign par_ini to either oldbestfit[[i]]$fitpars or if that doesn't exist, assign oldbestfit[[1]]$fitpars
  if (is.null(oldbestfit)) {
    par_ini_old = par_ini_full
  } else {
    if (i > length(oldbestfit)) {
      par_ini_old = oldbestfit[[1]]$fitpars
    } else {
      par_ini_old = oldbestfit[[i]]$fitpars
    }
  }

  replace_idx <- intersect(names(par_ini_full), names(par_ini_old))
  
  # replace starting values for fitted parameters with those from old best fit
  par_ini_full[replace_idx] <- par_ini_old[replace_idx]
  

  
  par_ini <- as.numeric(par_ini_full)

  # if we fit in log space, we transform
  # will be back-tronsformed inside fitting function before use in ODE model
  if (logfit==1) {
    par_ini <- log(par_ini)
    lb <- log(lb)
    ub <- log(ub)
  }

  # ---- fit ----
  bestfit <- nloptr::nloptr(
    x0 = par_ini,
    eval_f = fit_model1_function,
    lb = lb,
    ub = ub,
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
    tols = tols,
    simulatorname = simulatorname,
    logfit = logfit
  )
  # finished with fitting
  # doing some after fitting stuff

  # extract params
  params <- bestfit$solution
  names(params) <- fitparnames

# if we fit in log space, we transform back to get original values
  if (logfit==1) {
    params <- exp(params)
  }


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

# put the old and new objective function values together in a two-column data frame
old_objectives <- if (exists("oldbestfit")) {
  vapply(
    seq_along(bestfit_all),
    function(i) {
      if (!is.null(oldbestfit[[i]])) {
        oldbestfit[[i]]$objective
      } else {
        NA_real_
      }
    },
    numeric(1)
  )
} else {
  rep(NA_real_, length(bestfit_all))
}

new_objectives <- vapply(bestfit_all, function(x) x$objective, numeric(1))
objective_summary <- data.frame(
  old_objective = old_objectives,
  new_objective = new_objectives,
  improvement = old_objectives - new_objectives
)

print(objective_summary)


# play a sound when done
#beepr::beep(2)

# copy prior best fit to a new file if it exists
if (file.exists(here::here('results', 'output', bestfitfile))) {
  file.copy(
    from = here::here('results', 'output', bestfitfile),
    to = here::here('results', 'output', paste0('old',bestfitfile)),
    overwrite = TRUE
  )
  # save new best fit if it's better
  if (objective_summary$improvement[1] > 0) {
    saveRDS(bestfit_all, here::here('results', 'output', bestfitfile))
  }
} else {
  # no prior best fit, just save
  saveRDS(bestfit_all, here::here('results', 'output', bestfitfile))
}


