# -----------------------------------------------------------------------------
# model-config-function.R
# -----------------------------------------------------------------------------
# PURPOSE
#   Centralize model-specific configuration so the rest of the pipeline can
#   operate on a unified interface for either model1 or model2.
#
# WHAT THIS FILE CONTAINS
#   - Initial conditions (Y0)
#   - Initial parameter values and bounds
#   - Parameter labels for table output
#   - File paths for fixed parameters
#   - Simulator and objective function handles
# -----------------------------------------------------------------------------

library(here) # Project-root-relative file paths.

# Centralized virus transform helpers.
source(here::here("code", "analysis-code", "functions", "virus-transform-function.R"))

#' Build a model configuration object for model1 or model2.
#'
#' @param model_choice "model1" or "model2".
#' @return A list containing initial conditions, parameter defaults/bounds,
#'   labels, and the relevant simulator + fit functions.
#'
build_model_config <- function(model_choice) {
  # Validate input early to avoid confusing downstream errors.
  if (!model_choice %in% c("model1", "model2")) {
    stop("model_choice must be 'model1' or 'model2'.")
  }

  # Base parameter values and bounds (shared across models)
  par_ini_full <- c(
    b = 1e-8,
    k = 1e-5,
    p = 1e4,
    kF = 0.1,
    cV = 100,
    gF = 1,
    hV = 1e3,
    Fmax = 2,
    hF = 1,
    gS = 10,
    cS = 1,
    Emax_F = 1,
    C50_F = 1e-5,
    C50_V = 1e-8
  )

  lb <- c(
    b = 1e-12,
    k = 1e-10,
    p = 1,
    kF = 1e-1,
    cV = 0.1,
    gF = 1e-3,
    hV = 1e-2,
    Fmax = 0.1,
    hF = 1e-5,
    gS = 1e-3,
    cS = 1e-3,
    Emax_F = 1e-3,
    C50_F = 1e-10,
    C50_V = 1e-12
  )

  ub <- c(
    b = 1e-5,
    k = 1,
    p = 1e7,
    kF = 1e2,
    cV = 1e5,
    gF = 1e3,
    hV = 1e5,
    Fmax = 1e4,
    hF = 1e3,
    gS = 1e3,
    cS = 1e3,
    Emax_F = 1,
    C50_F = 1e-2,
    C50_V = 1e-2
  )

  # Descriptive labels used in tables/plots.
  parlabels_full <- c(
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
    sigma_add_IL6 = "Simga of IL6",
    sigma_add_WeightLossPerc = "Sigma of WeightLossPerc",
    sigma_prop_IL6 = "Proportional sigma of IL6",
    sigma_prop_WeightLossPerc = "Proportional sigma of WeightLossPerc"
  )

  # Add virus-specific sigma labels using the central quantity name.
  parlabels_full[paste0("sigma_add_", virus_quantity_name)] <- "Sigma of VirusLoad"
  parlabels_full[paste0("sigma_prop_", virus_quantity_name)] <- "Proportional sigma of VirusLoad"

  if (model_choice == "model1") {
    # Initial conditions for model 1 (no E compartment).
    Y0 <- c(
      Ad = 0,
      Ac = 0,
      At = 0,
      U = 1e7,
      I = 0,
      V = 1,
      F = 0,
      A = 0,
      S = 0
    )

    # Model 1 fixed parameters + simulator.
    fixedpars_file <- here::here("data", "processed-data", "model1-fixed-parameters.csv")
    simulatorname <- model1_simulator
    fit_function <- fit_model_function
  } else {
    # Initial conditions for model 2 (includes E compartment).
    Y0 <- c(
      Ad = 0,
      Ac = 0,
      At = 0,
      U = 1e7,
      E = 0,
      I = 0,
      V = 1,
      F = 0,
      A = 0,
      S = 0
    )

    # Model 2 fixed parameters + simulator.
    fixedpars_file <- here::here("data", "processed-data", "model2-fixed-parameters.csv")
    simulatorname <- model2_simulator
    fit_function <- fit_model_function
  }

  # Return a single configuration object used by all pipelines.
  list(
    Y0 = Y0,
    par_ini_full = par_ini_full,
    lb = lb,
    ub = ub,
    parlabels_full = parlabels_full,
    fixedpars_file = fixedpars_file,
    simulatorname = simulatorname,
    fit_function = fit_function
  )
}
