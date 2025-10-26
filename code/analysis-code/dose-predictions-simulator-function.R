#############################################################################
# function that runs the model over a wide range of dosing values
# and dosing schedules, and for each computes reduction
# in virus load, innate symptoms and weight loss
# it needs: deSolve (since it calls the simulator function), dplyr, caTools
#############################################################################
library(deSolve)
library(dplyr)
library(caTools)

# load model simulator function
source(here("code/analysis-code/model-simulator-function.R"))

simulate_dose_predictions <- function(bestfit) {
  # doses for which time-series is saved
  ts_doses <- c(1, 10, 100, 1000)

  # Model run parameters -------------------------------------------------------
  tfinal <- 7
  dt <- 0.005
  doses <- sort(unique(c(0, ts_doses, 10^seq(-3, 5, length = 100)))) #making sure we include ts_doses

  params <- bestfit$solution
  names(params) <- bestfit$fitparnames
  fixedpars <- bestfit$fixedpars
  Y0 <- bestfit$Y0

  #############################################################################
  # 2.  HELPER FUNCTIONS ------------------------------------------------------
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

  # Main simulator -------------------------------------------------------------
  simulate_dose_response <- function(
    doses,
    txstart,
    txend,
    txinterval,
    schedule_name,
    ts_doses
  ) {
    # initialize summary data frame
    summary_df <- tibble(
      Dose = doses,
      AUCV = NA_real_,
      AUCF = NA_real_,
      AUCS = NA_real_,
      Schedule = schedule_name
    )

    ts_store <- list() # collect full trajectories

    for (i in seq_along(doses)) {
      pars <- c(
        Y0,
        params,
        fixedpars,
        Ad0 = doses[i] / 50, # divide by body weight scaling to get actual dose
        txstart = txstart,
        txinterval = txinterval,
        txend = txend,
        tstart = 0,
        tfinal = tfinal,
        dt = dt
      )

      odeout <- try(do.call(simulate_model, as.list(pars)))
      if (inherits(odeout, "try-error") || length(odeout) == 1) {
        message("Integrator error – ", schedule_name, " (dose ", doses[i], ")")
        next
      }

      ode_df <- as.data.frame(odeout)

      # ---- AUCs --------------------------------------------------------------
      summary_df$AUCV[i] <- caTools::trapz(
        ode_df$time,
        log10(pmax(1, ode_df$V))
      )
      summary_df$AUCF[i] <- caTools::trapz(ode_df$time, ode_df$F)
      summary_df$AUCS[i] <- caTools::trapz(ode_df$time, ode_df$S)

      traj_df <- dplyr::mutate(
        ode_df,
        Dose = doses[i],
        Schedule = schedule_name
      )

      # ---- store time-series only for selected doses -------------------------
      if (doses[i] %in% ts_doses) {
        ts_store[[length(ts_store) + 1]] <- traj_df
      }
    }

    ret_list <- list(
      summary = summary_df,
      timeseries = if (length(ts_store)) bind_rows(ts_store) else tibble()
    )

    return(ret_list)
  }

  #############################################################################
  # 3.  RUN ALL SCHEDULES ------------------------------------------------------
  # each schedule is a different treatment regimen shown in the main text
  #############################################################################
  schedule_defs <- list(
    s1 = list(txstart = 1, txend = 4, txinterval = 0.5, name = "s1"), # baseline, as done in experiment
    s2 = list(txstart = 2, txend = 5, txinterval = 0.5, name = "s2"), # treatment start at day 2
    s3 = list(txstart = 3, txend = 6, txinterval = 0.5, name = "s3"), # treatment start at day 3
    s4 = list(txstart = 1, txend = 4, txinterval = 1, name = "s4"), # daily dosing
    s5 = list(txstart = 1, txend = 1, txinterval = 1, name = "s5") # single dosing
  )

  all_results <- lapply(
    schedule_defs,
    \(s) {
      simulate_dose_response(
        doses,
        s$txstart,
        s$txend,
        s$txinterval,
        s$name,
        ts_doses
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
    timeseries_df = timeseries_df
  )

  return(allres)
}
