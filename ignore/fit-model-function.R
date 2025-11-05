###################################################################
#function that fits the ODE model to data
#this function specifies the objective function that the nloptr solver is minimizing
#we are doing least squares fit here, but any function can be specified
#NLOPTR will try to minimize whatever is being fit
# it needs: nloptr, dplyr, deSolve (since it calls the simulator function)
library(nloptr)
library(dplyr)
library(deSolve)


###################################################################
# this function is set up to allow some parameters to be fitted (stored in params)
# and some parameters to be fixed (stored in fixedpars)
fit_model_function <- function(
  params,
  fitdata,
  Y0,
  tfinal,
  dt,
  fitparnames,
  fixedpars
) {
  #for some reason nloptr strips names from parameters
  #if i want to address parameters by name, I need to reassign their names
  names(params) = fitparnames

  Objfct_all = 0 #this will contain the objective function to be minimized

  # will contain all model fits
  allodeout = list()

  # compute maximum values for each variable, used to normalize when computing objective function
  #extract values for virus load at time points where data is available
  Vvals_all <- fitdata %>% filter(Quantity == "LogVirusLoad")
  Innvals_all <- fitdata %>% filter(Quantity == "IL6")
  Symvals_all <- fitdata %>% filter(Quantity == "WeightLossPerc")
  Vmax = max(Vvals_all$Value)
  Innmax = max(Innvals_all$Value)
  Symmax = max(Symvals_all$Value)

  # loop over 3 treatment scenarios
  # we fit all 3 scenarios at the same time here
  for (i in 1:length(doses)) {
    # combine parameters together to be sent to the simulate_ode function
    # treatment start/end/intervals are same for all scenarios
    allpars = c(
      Y0,
      params,
      fixedpars,
      Ad0 = doses[i],
      txstart = 1,
      txinterval = 0.5,
      txend = 4,
      tstart = 0,
      tfinal = tfinal,
      dt = dt
    )

    # this calls the simulate_ode function with the indicated parameters
    # the extra try function catches errors
    odeout <- try(do.call(simulate_model, as.list(allpars)))

    #try command catches error from ode function.
    # If error occurs and things "break", we exit the whole optimizer routine with a high objective function value,
    # this high value indicates that 'things didn't work'
    if (length(odeout) == 1) {
      cat(
        '!!!!!!unresolvable integrator error - triggering early return from optimizer!!!!!!'
      )
      return(1e10)
    }

    # save full fit object for later plotting
    allodeout[[i]] = odeout

    #extract values for virus load at time points where data is available
    Vvals <- Vvals_all %>% filter(Scenario == scenarios[i])
    Innvals <- Innvals_all %>% filter(Scenario == scenarios[i])
    Symvals <- Symvals_all %>% filter(Scenario == scenarios[i])

    # get model predictions for virus, log-transform, make sure no values below 0 (after transform)
    # any virus prediction that's below 1 virion will be set to 1 virion (so it's 0 after log transform)
    # that basically means any model predictions below 1 virion are treated as log(1)=0, same as in data
    # this means any time the model predicts values below LOD/0, it's set to that value
    Vpred = log10(pmax(1, odeout[match(Vvals$xvals, odeout[, "time"]), "V"]))
    #model predictions for innate response
    # both data and model predictions are in linear scale
    Innpred = odeout[match(Innvals$xvals, odeout[, "time"]), "F"]
    #model predictions for symptoms
    # both data and model predictions are in linear scale
    Sympred = odeout[match(Symvals$xvals, odeout[, "time"]), "S"]

    # compute sum of squares for each variable
    # weighted and normalized
    # this standardizes by the maximum of each variable across all scenarios/doses
    # standardizing for each scenario separately doesn't make sense, we just want to do it across variables

    Virtot = sum(((Vpred - Vvals$Value) / Vmax)^2) / nrow(Vvals)
    Inntot = sum(((Innpred - Innvals$Value) / Innmax)^2) / nrow(Innvals)
    Symtot = sum(((Sympred - Symvals$Value) / Symmax)^2) / nrow(Symvals)

    # objective function for a given treatment scenario
    # this is the sum of squares for each variable that is being fit
    Objfct_dose = Virtot + Inntot + Symtot

    # for diagnostics, if something is wrong with the objective function, we can inspect
    #if (is.nan(Objfct_dose)) {browser()}

    # compute overall objective function for all 3 treatment scenarios at the same time
    Objfct_all = Objfct_all + Objfct_dose
  } #end loop over doses

  return(Objfct_all) # return quantity to be minimized
} #end function that fits the ODE model to the data
