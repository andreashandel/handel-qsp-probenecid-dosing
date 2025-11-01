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
  scenarios
) {
  #for some reason nloptr strips names from parameters
  #if i want to address parameters by name, I need to reassign their names
  names(params) = fitparnames

  # --- Split parameters: structural vs sigma -------------------------------
  # Everything matching ^sigma_(add|prop)_* is treated as a sigma param.
  is_sigma <- grepl("^sigma_(add|prop)_", names(params))
  params_sigma <- params[is_sigma]
  params_struct <- params[!is_sigma]

  # For the ODE simulator, combine structural + fixed pars only
  ode_param_pool <- c(fixedpars, params_struct)

  # ---- Variance helpers (sigmas come ONLY from params_sigma) ----------------
  get_sigma <- function(qty, which = c("add", "prop")) {
    which <- match.arg(which)
    nm <- if (which == "add") {
      paste0("sigma_add_", qty)
    } else {
      paste0("sigma_prop_", qty)
    }
    if (!(nm %in% names(params_sigma))) {
      return(0)
    }
    val <- as.numeric(params_sigma[[nm]])
    if (!is.finite(val) || is.na(val)) {
      return(0)
    }
    abs(val)
  }

  # Pointwise variance: V = a^2 + (b * pred)^2
  var_fun <- function(pred, qty) {
    a <- get_sigma(qty, "add")
    b <- get_sigma(qty, "prop")
    v <- a^2 + (b * pmax(0, pred))^2
    pmax(v, 1e-12)
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
    allpars = c(
      Y0,
      ode_param_pool,
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
    #try command catches error from ode function.
    # If error occurs and things "break", we exit the whole optimizer routine with a high objective function value,
    # this high value indicates that 'things didn't work'
    odeout <- try(do.call(simulate_model, as.list(allpars)), silent = TRUE)
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
      if (nrow(Vvals)) {
        data.frame(
          Quantity = "LogVirusLoad",
          time = Vvals$xvals,
          obs = Vvals$Value,
          pred = Vpred
        )
      } else {
        NULL
      },
      if (nrow(Innvals)) {
        data.frame(
          Quantity = "IL6",
          time = Innvals$xvals,
          obs = Innvals$Value,
          pred = Innpred
        )
      } else {
        NULL
      },
      if (nrow(Symvals)) {
        data.frame(
          Quantity = "WeightLossPerc",
          time = Symvals$xvals,
          obs = Symvals$Value,
          pred = Sympred
        )
      } else {
        NULL
      }
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
      if (!nrow(di)) {
        next
      }

      Vi <- var_fun(di$pred, as.character(qty))
      res <- di$obs - di$pred

      nll_k <- 0.5 * sum(log(Vi) + (res^2) / Vi, na.rm = TRUE)
      w_k <- if (!is.null(w_by_qty[[as.character(qty)]])) {
        as.numeric(w_by_qty[[as.character(qty)]])
      } else {
        1
      }

      Objfct_dose <- Objfct_dose + w_k * nll_k
    }

    Objfct_all <- Objfct_all + Objfct_dose
  } # end loop over doses

  #browser()

  return(Objfct_all) # return quantity to be minimized
} #end function that fits the ODE model to the data
