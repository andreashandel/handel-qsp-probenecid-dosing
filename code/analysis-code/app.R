# Shiny application for interactively exploring the QSP model fits
#
# The app exposes the model parameters as numeric inputs so that users can
# explore how different parameter combinations affect the simulated time-series
# and the corresponding objective function value. The layout mirrors the static
# "best fit" figure by re-using the existing plotting helper.

# run with shiny::runApp("code/analysis-code/app.R")

library(shiny)
library(here)
library(dplyr)

# The plotting helper depends on ggplot2 and patchwork; loading them explicitly
# here ensures the app can render the figures even if the helper has not been
# sourced before.
library(ggplot2)
library(patchwork)

# Source the core model functions ------------------------------------------------
source(here::here("code", "analysis-code", "model-simulator-function.R"))
simulate_model_v1 <- simulate_model

source(here::here("code", "analysis-code", "model-simulator-function-v2.R"))
simulate_model_v2 <- simulate_model

# default binding keeps backwards compatibility for other scripts
simulate_model <- simulate_model_v1

source(here::here("code", "analysis-code", "fit-model-function.R"))
source(here::here("code", "plotting-code", "timeseries-plot-function.R"))

# Load and prepare the fitting data --------------------------------------------
fitdata_path <- here::here("data", "processed-data", "processeddata.csv")
fitdata <- read.csv(fitdata_path, stringsAsFactors = FALSE)

scenario_levels <- c("NoTreatment", "PanCytoVir10mg", "PanCytoVir100mg")
quantity_levels <- c("LogVirusLoad", "IL6", "WeightLossPerc")

fitdata <- fitdata |>
  mutate(
    Scenario = factor(Scenario, levels = scenario_levels),
    Quantity = factor(Quantity, levels = quantity_levels),
    Dose = c(0, 10, 100)[as.numeric(Scenario)],
    xvals = Day
  )

doses <- sort(unique(fitdata$Dose))
scenarios <- levels(fitdata$Scenario)

format_dose_label <- function(dose) {
  if (isTRUE(all.equal(dose, 0))) {
    return("no drug")
  }
  paste0(format(dose, trim = TRUE, scientific = FALSE), " mg/kg")
}

dose_labels <- vapply(doses, format_dose_label, character(1))

# Load fixed parameters ---------------------------------------------------------
fixedpars_path <- here::here("data", "processed-data", "fixed-parameters.csv")
fixedpars_df <- read.csv(fixedpars_path, stringsAsFactors = FALSE)
fixedpars_df$parname <- trimws(fixedpars_df$parname)
fixedpars_df$parnamefull <- trimws(fixedpars_df$parnamefull)

fixedpars_base <- fixedpars_df$value
names(fixedpars_base) <- fixedpars_df$parname
fixedpars_base <- as.numeric(fixedpars_base)
names(fixedpars_base) <- fixedpars_df$parname

# Initial conditions for the ODE system ----------------------------------------
Y0_base <- c(
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

# Baseline fitted parameters and bounds ----------------------------------------
fit_param_defaults <- c(
  b = 1e-8,
  k = 1e-4,
  p = 2e3,
  kF = 1,
  cV = 10,
  gF = 0.1,
  hV = 1e4,
  Fmax = 5,
  hF = 1,
  gS = 10,
  cS = 1,
  Emax_F = 0.5,
  C50_F = 1,
  C50_V = 1
)

fit_param_bounds <- data.frame(
  name = names(fit_param_defaults),
  min = c(1e-12, 1e-8, 1e-3, 1e-10, 1e-2, 1e-3, 1e-5, 0.1, 1e-3, 1e-4, 1e-2, 1e-3, 1e-7, 1e-7),
  max = c(1e-2, 1e5, 1e10, 1e3, 1e5, 1e3, 1e8, 10, 1e5, 1e4, 1e4, 1, 1e6, 1e6),
  stringsAsFactors = FALSE
)

# Measurement noise parameters --------------------------------------------------
var_by_qty <- fitdata |>
  group_by(Quantity) |>
  summarize(var_value = var(Value, na.rm = TRUE), .groups = "drop")

sigma_all <- c(
  sigma_add_LogVirusLoad = sqrt(var_by_qty$var_value[var_by_qty$Quantity == "LogVirusLoad"]),
  sigma_prop_LogVirusLoad = 0,
  sigma_add_IL6 = sqrt(var_by_qty$var_value[var_by_qty$Quantity == "IL6"]),
  sigma_prop_IL6 = 0,
  sigma_add_WeightLossPerc = sqrt(var_by_qty$var_value[var_by_qty$Quantity == "WeightLossPerc"]),
  sigma_prop_WeightLossPerc = 0
)

sigma_to_fit <- c(
  "sigma_add_LogVirusLoad",
  "sigma_add_IL6",
  "sigma_add_WeightLossPerc"
)

fit_param_defaults <- c(fit_param_defaults, sigma_all[sigma_to_fit])

sigma_bounds <- data.frame(
  name = names(sigma_all),
  min = 0,
  max = 1e3,
  stringsAsFactors = FALSE
)

fit_param_bounds <- bind_rows(
  fit_param_bounds,
  sigma_bounds[sigma_bounds$name %in% sigma_to_fit, ]
)
fit_param_bounds <- fit_param_bounds[!duplicated(fit_param_bounds$name), ]
fit_param_bounds$min <- 0

# Fixed parameters (including the non-fitted sigmas) ----------------------------
fixedpars_defaults <- c(fixedpars_base, sigma_all[setdiff(names(sigma_all), sigma_to_fit)])

# Attempt to load best-fit values as the starting point ------------------------
fit_param_defaults <- fit_param_defaults[!duplicated(names(fit_param_defaults))]
fixedpars_defaults <- fixedpars_defaults[!duplicated(names(fixedpars_defaults))]

fit_param_defaults_base <- fit_param_defaults
fixedpars_defaults_base <- fixedpars_defaults

load_bestfit_defaults <- function(bestfit_path, fit_defaults, fixed_defaults, Y0_defaults) {
  bestfit_obj <- tryCatch(readRDS(bestfit_path), error = function(e) NULL)

  if (is.list(bestfit_obj) && length(bestfit_obj) > 0) {
    bestfit_first <- bestfit_obj[[1]]

    if (!is.null(bestfit_first$fitpars)) {
      overlap <- intersect(names(fit_defaults), names(bestfit_first$fitpars))
      fit_defaults[overlap] <- bestfit_first$fitpars[overlap]
    } else if (!is.null(bestfit_first$solution) && !is.null(bestfit_first$fitparnames)) {
      bestfit_named <- bestfit_first$solution
      names(bestfit_named) <- bestfit_first$fitparnames
      overlap <- intersect(names(fit_defaults), names(bestfit_named))
      fit_defaults[overlap] <- bestfit_named[overlap]
    }

    if (!is.null(bestfit_first$fixedpars)) {
      overlap <- intersect(names(fixed_defaults), names(bestfit_first$fixedpars))
      fixed_defaults[overlap] <- bestfit_first$fixedpars[overlap]
    }

    if (!is.null(bestfit_first$Y0)) {
      overlap <- intersect(names(Y0_defaults), names(bestfit_first$Y0))
      Y0_defaults[overlap] <- bestfit_first$Y0[overlap]
    }
  }

  list(
    fit = fit_defaults[!duplicated(names(fit_defaults))],
    fixed = fixed_defaults[!duplicated(names(fixed_defaults))],
    Y0 = Y0_defaults[!duplicated(names(Y0_defaults))]
  )
}

build_model_config <- function(id, label, simulator, bestfit_path, Y0_override = NULL) {
  fit_defaults <- fit_param_defaults_base
  fixed_defaults <- fixedpars_defaults_base
  Y0_defaults <- Y0_base

  if (!is.null(Y0_override)) {
    overlap <- intersect(names(Y0_defaults), names(Y0_override))
    Y0_defaults[overlap] <- Y0_override[overlap]
  }

  defaults <- load_bestfit_defaults(bestfit_path, fit_defaults, fixed_defaults, Y0_defaults)

  if (!is.null(Y0_override)) {
    overlap <- intersect(names(defaults$Y0), names(Y0_override))
    defaults$Y0[overlap] <- Y0_override[overlap]
  }

  list(
    id = id,
    label = label,
    simulator = simulator,
    fit_defaults = defaults$fit,
    fixed_defaults = defaults$fixed,
    Y0 = defaults$Y0,
    bestfit_path = bestfit_path
  )
}

model_configs <- list(
  v1 = build_model_config(
    id = "v1",
    label = "model-simulator-function.R (A = 0)",
    simulator = simulate_model_v1,
    bestfit_path = here::here("results", "output", "bestfit.Rds"),
    Y0_override = c(A = 0)
  ),
  v2 = build_model_config(
    id = "v2",
    label = "model-simulator-function-v2.R (A = 1)",
    simulator = simulate_model_v2,
    bestfit_path = here::here("results", "output", "bestfit-v2.Rds"),
    Y0_override = c(A = 1)
  )
)

default_model_id <- names(model_configs)[1]
initial_fit_defaults <- model_configs[[default_model_id]]$fit_defaults
initial_fixed_defaults <- model_configs[[default_model_id]]$fixed_defaults
fit_param_names <- names(fit_param_defaults_base)
fixed_param_names <- names(fixedpars_defaults_base)
model_choices <- setNames(
  names(model_configs),
  vapply(model_configs, function(cfg) cfg$label, character(1))
)

# Helper metadata for UI construction -----------------------------------------
fit_param_labels <- c(
  b = "Virus infection rate",
  k = "Adaptive response clearance rate",
  p = "Virus production rate",
  kF = "Innate response suppression strength",
  cV = "Virus removal rate",
  gF = "Maximum innate response induction",
  hV = "Half-saturation for virus-induced innate activation",
  Fmax = "Maximum innate response",
  hF = "Adaptive response half-maximum induction",
  gS = "Symptom induction rate",
  cS = "Symptom decay rate",
  Emax_F = "Maximum innate response suppression",
  C50_F = "Half maximum of innate response effect",
  C50_V = "Half maximum of virus suppression effect",
  sigma_add_LogVirusLoad = "Sigma (additive) – Log Virus Load",
  sigma_add_IL6 = "Sigma (additive) – IL-6",
  sigma_add_WeightLossPerc = "Sigma (additive) – Weight loss"
)

fixed_param_labels <- fixedpars_df$parnamefull
names(fixed_param_labels) <- fixedpars_df$parname

sigma_fixed_labels <- c(
  sigma_prop_LogVirusLoad = "Sigma (proportional) – Log Virus Load",
  sigma_prop_IL6 = "Sigma (proportional) – IL-6",
  sigma_prop_WeightLossPerc = "Sigma (proportional) – Weight loss"
)

# Re-usable helper to evaluate the model --------------------------------------
solvertype <- "vode"
tols <- 1e-9
tfinal <- 7
dt <- 0.02

param_input_id <- function(name) paste0("par_", name)
fixed_input_id <- function(name) paste0("fixed_", name)

`%||%` <- function(x, y) if (is.null(x)) y else x

format_input_label <- function(name, label_map) {
  label_text <- label_map[[name]]
  label_text <- if (is.null(label_text)) "" else as.character(label_text)

  if (!nzchar(label_text)) {
    return(name)
  }

  paste(name, label_text, sep = " - ")
}

format_condition <- function(cond) {
  if (inherits(cond, "condition")) {
    conditionMessage(cond)
  } else if (is.null(cond)) {
    ""
  } else {
    as.character(cond)
  }
}

run_model_once <- function(params, fixedpars, Y0_vals, simulator) {
  params <- params[!is.na(params)]
  fixedpars <- fixedpars[!is.na(fixedpars)]

  params_ode <- params[!grepl("^sigma_", names(params))]
  fixedpars_ode <- fixedpars[!grepl("^sigma_", names(fixedpars))]

  simulate_one <- function(ad0, scenario_label) {
    args <- c(
      as.list(Y0_vals),
      as.list(params_ode),
      as.list(fixedpars_ode),
      list(
        Ad0 = ad0,
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

    odeout <- tryCatch(do.call(simulator, args), error = identity)

    if (inherits(odeout, "error")) {
      return(odeout)
    }

    df <- as.data.frame(odeout)
    df$Dose <- ad0
    df$Scenario <- factor(scenario_label, levels = scenarios)
    df
  }

  sim_list <- lapply(seq_along(doses), function(i) {
    simulate_one(doses[i], scenarios[i])
  })

  failed <- vapply(sim_list, inherits, logical(1), "error")
  if (any(failed)) {
    first_error <- sim_list[[which(failed)[1]]]
    return(list(error = first_error))
  }

  sim_df <- bind_rows(sim_list)

  objective <- tryCatch(
    fit_model_function(
      params = unname(params),
      fitdata = fitdata,
      Y0 = Y0_vals,
      tfinal = tfinal,
      dt = dt,
      fitparnames = names(params),
      fixedpars = fixedpars,
      doses = doses,
      scenarios = scenarios,
      solvertype = solvertype,
      tols = tols,
      simulator = simulator
    ),
    error = identity
  )

  list(sim = sim_df, objective = objective)
}

# ---------------------------------------------------------------------------- #
# Shiny user interface ---------------------------------------------------------
# ---------------------------------------------------------------------------- #
ui <- fluidPage(
  titlePanel("PanCytoVir model explorer"),
  tags$head(
    tags$style(
      HTML(
        ".parameter-grid {\n           display: grid;\n           grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));\n           gap: 0.75rem;\n         }\n         .parameter-grid .shiny-input-container {\n           margin-bottom: 0;\n         }\n         .control-header {\n           margin-top: 1rem;\n           font-weight: 600;\n         }"
      )
    )
  ),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      selectInput(
        "model_choice",
        "Model version",
        choices = model_choices,
        selected = default_model_id
      ),
      helpText(
        "Adjust the parameters to immediately update the simulation ",
        "and objective function evaluation."
      ),
      br(),
      tabsetPanel(
        id = "param-tabs",
        tabPanel(
          "Fitted parameters",
          div(
            class = "parameter-grid",
            lapply(fit_param_names, function(nm) {
              bounds_row <- fit_param_bounds[fit_param_bounds$name == nm, ]
              max_val <- if (nrow(bounds_row)) bounds_row$max else NA
              label <- format_input_label(nm, fit_param_labels)
              numericInput(
                inputId = param_input_id(nm),
                label = label,
                value = initial_fit_defaults[[nm]],
                min = 0,
                max = max_val,
                step = NA
              )
            })
          )
        ),
        tabPanel(
          "Fixed parameters",
          div(
            class = "parameter-grid",
            lapply(names(fixedpars_base), function(nm) {
              label <- format_input_label(nm, fixed_param_labels)
              numericInput(
                inputId = fixed_input_id(nm),
                label = label,
                value = initial_fixed_defaults[[nm]],
                min = 0,
                step = NA
              )
            })
          ),
          div(
            class = "parameter-grid",
            tags$div(class = "control-header", "Measurement noise (fixed)"),
            lapply(setdiff(names(sigma_all), sigma_to_fit), function(nm) {
              label <- format_input_label(nm, sigma_fixed_labels)
              sigma_max <- sigma_bounds$max[sigma_bounds$name == nm]
              numericInput(
                inputId = fixed_input_id(nm),
                label = label,
                value = initial_fixed_defaults[[nm]],
                min = 0,
                max = ifelse(length(sigma_max), sigma_max[1], NA),
                step = NA
              )
            })
          )
        )
      )
    ),
    mainPanel(
      width = 8,
      h4("Objective function"),
      textOutput("objective_value"),
      br(),
      plotOutput("fit_plot", height = "900px")
    )
  )
)

# ---------------------------------------------------------------------------- #
# Server logic -----------------------------------------------------------------
# ---------------------------------------------------------------------------- #
server <- function(input, output, session) {
  model_config <- reactive({
    req(input$model_choice)
    model_configs[[input$model_choice]]
  })

  observeEvent(model_config(), {
    cfg <- model_config()

    for (nm in fit_param_names) {
      updateNumericInput(
        session = session,
        inputId = param_input_id(nm),
        value = cfg$fit_defaults[[nm]]
      )
    }

    for (nm in fixed_param_names) {
      updateNumericInput(
        session = session,
        inputId = fixed_input_id(nm),
        value = cfg$fixed_defaults[[nm]]
      )
    }
  }, ignoreNULL = FALSE)

  fitted_values <- reactive({
    cfg <- model_config()
    defaults <- cfg$fit_defaults

    vapply(
      fit_param_names,
      function(nm) input[[param_input_id(nm)]] %||% defaults[[nm]],
      numeric(1),
      USE.NAMES = TRUE
    )
  })

  fixed_values <- reactive({
    cfg <- model_config()
    defaults <- cfg$fixed_defaults

    vapply(
      fixed_param_names,
      function(nm) input[[fixed_input_id(nm)]] %||% defaults[[nm]],
      numeric(1),
      USE.NAMES = TRUE
    )
  })

  model_results <- reactive({
    cfg <- model_config()
    params <- fitted_values()
    fixedpars <- fixed_values()

    run_model_once(params, fixedpars, cfg$Y0, cfg$simulator)
  })

  output$objective_value <- renderText({
    res <- model_results()

    if (!is.null(res$error)) {
      return(paste("Simulation failed:", format_condition(res$error)))
    }

    if (inherits(res$objective, "error")) {
      return(paste("Objective evaluation failed:", format_condition(res$objective)))
    }

    paste0("Current objective value: ", format(res$objective, digits = 6, scientific = TRUE))
  })

  output$fit_plot <- renderPlot({
    res <- model_results()
    validate(need(is.null(res$error), format_condition(res$error)))
    validate(need(!inherits(res$objective, "error"), format_condition(res$objective)))

    sim_df <- res$sim |>
      mutate(
        Dose = factor(Dose, levels = doses)
      )

    plot_timeseries(
      data = fitdata,
      modelfit = sim_df,
      tmax = tfinal,
      dose_levels = doses,
      dose_levels_labels = dose_labels,
      x_jitter = 0.3
    )
  }, res = 96)
}

if (interactive()) {
  shinyApp(ui, server)
}
