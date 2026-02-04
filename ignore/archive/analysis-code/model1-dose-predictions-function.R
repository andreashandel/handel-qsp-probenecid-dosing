#############################################################################
# function that runs the fitted model over a wide range of dosing values
# and dosing schedules, and for each computes reduction
# in virus load, innate symptoms and weight loss
# it needs: deSolve (since it calls the simulator function), dplyr, caTools
#############################################################################
library(deSolve)
library(dplyr)
library(caTools)
library(here)

# load model simulator function
source(here("code/analysis-code/model1-simulator-function.R"))

simulate_dose_predictions <- function(bestfit, ts_doses, solvertype, tols, dt, tfinal, simulatorname) 
{
  
  
  #############################################################################
  #  Inner FUNCTIONS ------------------------------------------------------
  #############################################################################

  # Percent-reduction helper ---------------------------------------------------
  compute_percent_reduction <- function(df) {
    df %>%
      group_by(Schedule) %>%
      mutate(
        perc_AUCV = (first(AUCV) - AUCV) / first(AUCV) * 100,
        perc_AUCF = (first(AUCF) - AUCF) / first(AUCF) * 100,
        perc_AUCS = (first(AUCS) - AUCS) / first(AUCS) * 100
      )
  }


  # Main simulator function-------------------------------------------------------------
  simulate_dose_response <- function(
    all_doses,
    txstart,
    txend,
    txinterval,
    schedule_name,
    ts_doses,
    Y0,
    fitpars_ode,
    fixedpars_ode,
    solvertype,
    tols,
    dt, 
    tfinal
  ) {
    # initialize summary data frame
    summary_df <- tibble(
      Dose = all_doses,
      AUCV = NA_real_,
      AUCF = NA_real_,
      AUCS = NA_real_,
      Schedule = schedule_name
    )

    ts_store <- list() # collect full trajectories

    for (i in seq_along(all_doses)) {
    allpars <- c(
      as.list(Y0),
      as.list(fitpars_ode),
      as.list(fixedpars_ode),
      list(
        Ad0 = all_doses[i],
        txstart = txstart,
        txinterval = txinterval,
        txend = txend,
        tstart = 0,
        tfinal = tfinal,
        dt = dt,
        solvertype = solvertype,
        tols = tols
      )
    )
      
    # odeout <- try(do.call(simulate_model, allpars), silent = TRUE)
    #   if (inherits(odeout, "try-error") || length(odeout) == 1) {
    #     message("Integrator error – ", schedule_name, " (dose ", all_doses[i], ")")
    #     next
    #   }
    
    odeout <- do.call(simulatorname, allpars)
    


    ode_df <- as.data.frame(odeout)
      
    # add free drug and PD effect variables for plotting; mirror ODE definitions
    Ct <- ode_df$At / as.numeric(allpars["Vt"])
    fu <- as.numeric(allpars["fmax"]) * Ct / (as.numeric(allpars["f50"]) + Ct)
    Cu <- fu * Ct
    fV <- as.numeric(allpars["Emax_V"]) * Cu / (as.numeric(allpars["C50_V"]) + Cu) # effect of drug on virus
    fF <- as.numeric(allpars["Emax_F"]) * Cu / (as.numeric(allpars["C50_F"]) + Cu) # effect of drug on innate
    ode_df$Cu <- Cu 
    ode_df$fV <- fV
    ode_df$fF <- fF 
       

      # ---- AUCs --------------------------------------------------------------
      summary_df$AUCV[i] <- caTools::trapz(
        ode_df$time,
        log10(pmax(1, ode_df$V))
      )
      summary_df$AUCF[i] <- caTools::trapz(ode_df$time, ode_df$F)
      summary_df$AUCS[i] <- caTools::trapz(ode_df$time, ode_df$S)

      traj_df <- dplyr::mutate(
        ode_df,
        Dose = all_doses[i],
        Schedule = schedule_name
      )

      # ---- store time-series only for selected doses -------------------------
      if (all_doses[i] %in% ts_doses) {
        ts_store[[length(ts_store) + 1]] <- traj_df
      }
    }

    ret_list <- list(
      summary = summary_df,
      timeseries = if (length(ts_store)) bind_rows(ts_store) else tibble()
    )

    return(ret_list)
  }


  # Main function part, calls the functions above
  # Model run parameters -------------------------------------------------------
  tfinal <- 7
  all_doses <- sort(unique(c(ts_doses, 10^seq(-2, 5, length = 100)))) #making sure we include ts_doses

  #params <- bestfit$solution
  #names(params) <- bestfit$fitparnames
  fitpars <- bestfit$fitpars
  #names(params) <- bestfit$fitparnames
  fixedpars <- bestfit$fixedpars
  Y0 <- bestfit$Y0

# pull out fitted parameters that are part of the ODE, excluding the sigmas
  fit_sigmas <- grepl("^sigma_(add|prop)_", names(fitpars))
  fitpars_ode = fitpars[!fit_sigmas]

  # pull out fixed parameters that are part of the ODE, excluding the sigmas
  fixed_sigmas <- grepl("^sigma_(add|prop)_", names(fixedpars))
  fixedpars_ode = fixedpars[!fixed_sigmas]


  
  #############################################################################
  # 3.  RUN ALL SCHEDULES ------------------------------------------------------
  # each schedule is a different treatment regimen shown in the main text
  #############################################################################
  schedule_defs <- list(
    s1 = list(txstart = 1, txend = 3.9, txinterval = 0.5, name = "s1"), # baseline, as done in experiment
    s2 = list(txstart = 2, txend = 4.9, txinterval = 0.5, name = "s2"), # treatment start at day 2
    s3 = list(txstart = 3, txend = 5.9, txinterval = 0.5, name = "s3"), # treatment start at day 3
    s4 = list(txstart = 1, txend = 3.9, txinterval = 1, name = "s4"), # daily dosing
    s5 = list(txstart = 1, txend = 1, txinterval = 1, name = "s5") # single dosing
  )

  all_results <- lapply(
    schedule_defs,
    \(s) {
      simulate_dose_response(
        all_doses,
        s$txstart,
        s$txend,
        s$txinterval,
        s$name,
        ts_doses,
        Y0,
    fitpars_ode,
    fixedpars_ode,
    solvertype,
    tols,
    dt,
    tfinal

      )
    }
  )

  # Split summary vs. trajectories --------------------------------------------
  all_results_df <- bind_rows(lapply(all_results, `[[`, "summary"))
  timeseries_df <- bind_rows(lapply(all_results, `[[`, "timeseries"))

  #############################################################################
  # 4.  POST-PROCESS & RETURN ---------------------------------------------------
  #############################################################################
  label_map <- c(
    s1 = "baseline",
    s2 = "d2 start",
    s3 = "d3 start",
    s4 = "daily tx",
    s5 = "single tx"
  )

  reduction_df <- compute_percent_reduction(all_results_df) %>%
    mutate(Scenario = label_map[Schedule])

  allres = list(
    all_results_df = all_results_df,
    reduction_df = reduction_df,
    timeseries_df = timeseries_df,
    ts_doses = ts_doses
  )

  return(allres)
}
