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
  fixedpars,
  doses,
  scenarios,
  solvertype,
  tols,
  simulator = simulate_model
) {
  simulator <- match.fun(simulator)

  #for some reason nloptr strips names from parameters
  #if i want to address parameters by name, I need to reassign their names
  names(params) = fitparnames

  # pull out fitted parameters that are part of the ODE, excluding the sigmas
  fit_sigmas <- grepl("^sigma_(add|prop)_", names(params))
  fitpars_ode = params[!fit_sigmas]

  # pull out fixed parameters that are part of the ODE, excluding the sigmas
  fixed_sigmas <- grepl("^sigma_(add|prop)_", names(fixedpars))
  fixedpars_ode = fixedpars[!fixed_sigmas]

  # Build a sigma pool that includes both FITTED and FIXED sigma parameters
  sigma_pool <- c(params[fit_sigmas], fixedpars[fixed_sigmas])

  # Helpers to fetch sigma values from the combined pool
  get_sigma <- function(qty, which = c("add", "prop")) {
    which <- match.arg(which)
    nm <- paste0("sigma_", which, "_", qty)
    as.numeric(sigma_pool[[nm]])
  }

  # Pointwise variance: V = a^2 + (b * pred)^2
  var_fun <- function(pred, qty) {
    a <- get_sigma(qty, "add")
    b <- get_sigma(qty, "prop")
    v <- a^2 + (b * pred)^2
  }

  # ---- Per-quantity weights to equalize contributions -----------------------
  n_by_qty <- fitdata %>%
    group_by(Quantity) %>%
    summarize(n = n(), .groups = "drop")
  w_by_qty <- setNames(1 / n_by_qty$n, n_by_qty$Quantity)
  w_by_qty <- w_by_qty / mean(w_by_qty, na.rm = TRUE)

  Objfct_all <- 0

  # will contain all model fits
  #allodeout <- vector("list", length(doses))

  # loop over all treatment scenarios
  # we fit all scenarios at the same time here
  for (i in seq_along(doses)) {
    # combine parameters together to be sent to the simulate_ode function
    # treatment start/end/intervals are same for all scenarios
    allpars <- c(
      as.list(Y0),
      as.list(fitpars_ode),
      as.list(fixedpars_ode),
      list(
        Ad0 = doses[i],
        txstart = 1,
        txinterval = 0.5,
        txend = 4,
        tstart = 0,
        tfinal = tfinal,
        dt = dt,
        solvertype = solvertype,
        tols = tols
      )
    )

    # this calls the simulate_ode function with the indicated parameters
    # the extra try function catches errors
    #try command catches error from ode function.
    # If error occurs and things "break", we exit the whole optimizer routine with a high objective function value,
    # this high value indicates that 'things didn't work'
    odeout <- try(do.call(simulator, allpars), silent = TRUE)
    if (inherits(odeout, "try-error")) {
      cat(
        "!!!!!!unresolvable integrator error - triggering early return from optimizer!!!!!!"
      )
      return(1e10)
    }

    # save full fit object for later plotting
    # allodeout[[i]] = odeout

    #extract values for virus load at time points where data is available
    Vvals <- fitdata %>%
      filter(Quantity == "LogVirusLoad", Scenario == scenarios[i])
    Innvals <- fitdata %>% filter(Quantity == "IL6", Scenario == scenarios[i])
    Symvals <- fitdata %>%
      filter(Quantity == "WeightLossPerc", Scenario == scenarios[i])

    tvec <- odeout[, "time"]

    # get model predictions for virus, log-transform, make sure no values below 0 (after transform)
    # any virus prediction that's below 1 virion will be set to 1 virion (so it's 0 after log transform)
    # that basically means any model predictions below 1 virion are treated as log(1)=0, same as in data
    # this means any time the model predicts values below LOD/0, it's set to that value
    # predictions, aligned to observed times (with interpolation)
    Vpred = log10(pmax(1, odeout[match(Vvals$xvals, tvec), "V"]))
    #model predictions for innate response
    # both data and model predictions are in linear scale
    Innpred = odeout[match(Innvals$xvals, tvec), "F"]
    #model predictions for symptoms
    # both data and model predictions are in linear scale
    Sympred = odeout[match(Symvals$xvals, tvec), "S"]

    pieces <- list(
      data.frame(
        Quantity = "LogVirusLoad",
        time = Vvals$xvals,
        obs = Vvals$Value,
        pred = Vpred
      ),
      data.frame(
        Quantity = "IL6",
        time = Innvals$xvals,
        obs = Innvals$Value,
        pred = Innpred
      ),
      data.frame(
        Quantity = "WeightLossPerc",
        time = Symvals$xvals,
        obs = Symvals$Value,
        pred = Sympred
      )
    )
    pieces <- Filter(NROW, pieces)
    if (!length(pieces)) {
      next
    }
    df <- do.call(rbind, pieces)

    # Contribution by variable, using chosen error model and sigmas
    Objfct_dose <- 0
    for (qty in unique(df$Quantity)) {
      di <- df[df$Quantity == qty, , drop = FALSE]

      Vi <- var_fun(di$pred, as.character(qty))
      res <- di$obs - di$pred

      nll_k <- 0.5 * sum(log(Vi) + (res^2) / Vi)
      w_k <- as.numeric(w_by_qty[[as.character(qty)]])

      Objfct_dose <- Objfct_dose + w_k * nll_k
    }

    Objfct_all <- Objfct_all + Objfct_dose
  } # end loop over doses

  return(Objfct_all) # return quantity to be minimized
} #end function that fits the ODE model to the data
