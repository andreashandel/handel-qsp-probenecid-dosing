# -----------------------------------------------------------------------------
# app.R
# -----------------------------------------------------------------------------
# PURPOSE
#   Shiny application for interactively exploring the QSP model fits.
#   Users can adjust parameters, view simulated time-series, and track the
#   objective function value in real time.
#
# HIGH-LEVEL FLOW (READ ME FIRST)
#   1) Load data, fixed parameters, and sigma defaults used by the fitter.
#   2) Build model configs (model1/model2) with defaults + bestfit overrides.
#   3) Construct the UI (parameter inputs + plot output).
#   4) On every input change, simulate the model and compute the objective.
#   5) Render the time-series plot and objective value for the user.
# -----------------------------------------------------------------------------
#
# DETAILED WALKTHROUGH
#   - Data setup: uses load_fit_data() for consistent factor levels and doses.
#   - Parameter defaults: pulls baseline values from model-config and fixed CSVs.
#   - Bestfit overrides: uses previously saved bestfit RDS if present.
#   - UI: two tabs for fitted vs fixed parameters, plus a plot panel.
#   - Server: reads inputs, simulates trajectories, and computes the objective.

library(shiny)     # Web app framework.
library(here)      # Project-root-relative file paths.
library(dplyr)     # Data manipulation utilities.
library(tidyr)     # Pivoting for objective contribution table.
library(ggplot2)   # Base plotting for time-series output.
library(patchwork) # Plot layout helper used in timeseries plotting.

# Centralized virus transform helpers.
source(here::here("code", "analysis-code", "functions", "virus-transform-function.R"))

# Model simulators + fit function
source(here::here("code", "analysis-code", "functions", "model1-simulator-function.R"))
source(here::here("code", "analysis-code", "functions", "model2-simulator-function.R"))
source(here::here("code", "analysis-code", "functions", "fit-function.R"))

# Shared helpers
source(here::here("code", "analysis-code", "functions", "fit-data-function.R"))
source(here::here("code", "analysis-code", "functions", "fixed-parameters-function.R"))
source(here::here("code", "analysis-code", "functions", "sigma-settings-function.R"))
source(here::here("code", "analysis-code", "functions", "model-config-function.R"))
source(here::here("code", "analysis-code", "functions", "objective-components-function.R"))

# Plotting helper (for time-series layout)
source(here::here("code", "plotting-code", "functions", "timeseries-plot-function.R"))

# -----------------------------------------------------------------------------
# Load and prepare data
# -----------------------------------------------------------------------------
fitdata_bundle <- load_fit_data()
fitdata <- fitdata_bundle$fitdata
scenarios <- fitdata_bundle$scenarios
doses <- fitdata_bundle$doses

format_dose_label <- function(dose) {
  if (isTRUE(all.equal(dose, 0))) {
    return("no drug")
  }
  paste0(format(dose, trim = TRUE, scientific = FALSE), " mg/kg")
}

dose_labels <- vapply(doses, format_dose_label, character(1))

# -----------------------------------------------------------------------------
# Fixed parameters + sigma settings (shared across models)
# -----------------------------------------------------------------------------
fixedpars_model1 <- load_fixed_parameters(
  here::here("data", "processed-data", "model1-fixed-parameters.csv")
)
fixedpars_model2 <- load_fixed_parameters(
  here::here("data", "processed-data", "model2-fixed-parameters.csv")
)

weight_mode <- "equal_quantity"

objective_data <- prepare_fast_objective_data(fitdata, weight_mode = weight_mode)
obs_times_by_scenario <- lapply(objective_data$data_by_scenario, function(x) {
  sort(unique(c(0, x$time)))
})

sigma_settings <- compute_sigma_settings(
  fitdata,
  sigma_to_fit = character(0)
)

sigma_fixed_labels <- c(
  setNames("Sigma (proportional) – Virus Load", paste0("sigma_prop_", virus_quantity_name)),
  sigma_prop_IL6 = "Sigma (proportional) – IL-6",
  sigma_prop_WeightLossPerc = "Sigma (proportional) – Weight loss"
)

# -----------------------------------------------------------------------------
# Base model configs (shared defaults, updated with bestfit if available)
# -----------------------------------------------------------------------------
select_bestfit_entry <- function(bestfit_obj) {
  if (!is.list(bestfit_obj) || !length(bestfit_obj)) {
    return(NULL)
  }

  # If this already looks like a packed bestfit, return it directly.
  if (!is.null(bestfit_obj$fitpars) || !is.null(bestfit_obj$solution)) {
    return(bestfit_obj)
  }

  objectives <- vapply(bestfit_obj, function(entry) {
    if (is.list(entry) && !is.null(entry$objective) && is.finite(entry$objective)) {
      entry$objective
    } else {
      Inf
    }
  }, numeric(1))

  if (all(!is.finite(objectives))) {
    return(bestfit_obj[[1]])
  }

  bestfit_obj[[which.min(objectives)]]
}

resolve_bestfit_path <- function(bestfit_paths) {
  for (path in bestfit_paths) {
    if (is.character(path) && length(path) == 1 && file.exists(path)) {
      return(path)
    }
  }
  NULL
}

load_bestfit_defaults <- function(bestfit_paths, fit_defaults, fixed_defaults, Y0_defaults) {
  bestfit_path <- resolve_bestfit_path(bestfit_paths)
  if (is.null(bestfit_path)) {
    return(list(
      fit = fit_defaults[!duplicated(names(fit_defaults))],
      fixed = fixed_defaults[!duplicated(names(fixed_defaults))],
      Y0 = Y0_defaults[!duplicated(names(Y0_defaults))]
    ))
  }

  bestfit_obj <- tryCatch(readRDS(bestfit_path), error = function(e) NULL)

  if (is.list(bestfit_obj) && length(bestfit_obj) > 0) {
    bestfit_first <- select_bestfit_entry(bestfit_obj)
    if (is.null(bestfit_first)) {
      bestfit_first <- bestfit_obj[[1]]
    }

    if (!is.null(bestfit_first$fitpars)) {
      extra <- setdiff(names(bestfit_first$fitpars), names(fit_defaults))
      if (length(extra)) {
        fit_defaults <- c(fit_defaults, bestfit_first$fitpars[extra])
      }
      overlap <- intersect(names(fit_defaults), names(bestfit_first$fitpars))
      fit_defaults[overlap] <- bestfit_first$fitpars[overlap]
    } else if (!is.null(bestfit_first$solution) && !is.null(bestfit_first$fitparnames)) {
      bestfit_named <- bestfit_first$solution
      names(bestfit_named) <- bestfit_first$fitparnames
      extra <- setdiff(names(bestfit_named), names(fit_defaults))
      if (length(extra)) {
        fit_defaults <- c(fit_defaults, bestfit_named[extra])
      }
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

build_app_model_config <- function(model_choice, fixedpars_bundle, bestfit_path) {
  config <- build_model_config(model_choice)

  fit_defaults <- c(config$par_ini_full, sigma_settings$sigma_fit_ini)
  fixed_defaults <- c(fixedpars_bundle$values, sigma_settings$sigma_fixed)

  defaults <- load_bestfit_defaults(bestfit_path, fit_defaults, fixed_defaults, config$Y0)

  list(
    id = model_choice,
    label = ifelse(model_choice == "model1", "Model 1", "Model 2"),
    simulatorname = config$simulatorname,
    fit_defaults = defaults$fit,
    fixed_defaults = defaults$fixed,
    Y0 = defaults$Y0,
    bestfit_path = bestfit_path,
    fit_param_names = names(defaults$fit),
    fixed_param_names = names(fixedpars_bundle$values),
    fixed_param_labels = fixedpars_bundle$labels
  )
}

model_configs <- list(
  model1 = build_app_model_config(
    model_choice = "model1",
    fixedpars_bundle = fixedpars_model1,
    bestfit_path = c(
      here::here("results", "output", "model1-bestfit-multistart.Rds"),
      here::here("results", "output", "model1-bestfit-sample.Rds")
    )
  ),
  model2 = build_app_model_config(
    model_choice = "model2",
    fixedpars_bundle = fixedpars_model2,
    bestfit_path = c(
      here::here("results", "output", "model2-bestfit-multistart.Rds"),
      here::here("results", "output", "model2-bestfit-sample.Rds")
    )
  )
)

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
  setNames("Sigma (additive) – Virus Load", paste0("sigma_add_", virus_quantity_name)),
  sigma_add_IL6 = "Sigma (additive) – IL-6",
  sigma_add_WeightLossPerc = "Sigma (additive) – Weight loss"
)

# -----------------------------------------------------------------------------
# UI helpers
# -----------------------------------------------------------------------------
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

format_breakdown_table <- function(breakdown_df, value_col) {
  if (is.null(breakdown_df)) {
    return(NULL)
  }

  if (is.list(breakdown_df) && !is.data.frame(breakdown_df)) {
    breakdown_df <- dplyr::bind_rows(breakdown_df)
  }

  if (!is.data.frame(breakdown_df)) {
    return(NULL)
  }

  required_cols <- c("Quantity", "Scenario", value_col)
  if (!all(required_cols %in% names(breakdown_df))) {
    return(NULL)
  }

  breakdown_df %>%
    mutate(
      Quantity = dplyr::recode(
        as.character(Quantity),
        VirusLoad = "V",
        IL6 = "F",
        WeightLossPerc = "S",
        .default = as.character(Quantity)
      ),
      Scenario = as.character(Scenario)
    ) %>%
    group_by(Quantity, Scenario) %>%
    summarise(value = sum(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = Scenario,
      values_from = value,
      values_fill = NA_real_,
      values_fn = sum
    )
}

compute_fast_weighted_nll_breakdown <- function(sim_obj_df, sigma_pool, objective_data) {
  quantity_levels <- objective_data$quantity_levels
  quantity_id_map <- objective_data$quantity_id_map
  scenario_levels <- objective_data$scenario_levels

  sigma_add_by_id <- rep(0, length(quantity_levels))
  sigma_prop_by_id <- rep(0, length(quantity_levels))

  for (i in seq_along(quantity_levels)) {
    qname <- quantity_levels[i]
    add_name <- paste0("sigma_add_", qname)
    prop_name <- paste0("sigma_prop_", qname)
    if (add_name %in% names(sigma_pool)) {
      sigma_add_by_id[i] <- sigma_pool[[add_name]]
    }
    if (prop_name %in% names(sigma_pool)) {
      sigma_prop_by_id[i] <- sigma_pool[[prop_name]]
    }
  }

  virus_id <- quantity_id_map[[virus_quantity_name]]
  il6_id <- quantity_id_map[["IL6"]]
  wl_id <- quantity_id_map[["WeightLossPerc"]]

  block_rows <- lapply(names(objective_data$data_by_scenario), function(scenario_key) {
    scenario_data <- objective_data$data_by_scenario[[scenario_key]]
    ode_df <- sim_obj_df %>%
      filter(as.character(Scenario) == scenario_key) %>%
      select(time, V, F, S)

    time_index <- match(scenario_data$time, ode_df$time)
    if (any(is.na(time_index))) {
      return(NULL)
    }

    qid <- scenario_data$quantity_id
    pred <- numeric(length(qid))

    if (!is.null(virus_id) && any(qid == virus_id)) {
      v_pred <- transform_virus(ode_df$V[time_index])
      pred[qid == virus_id] <- v_pred[qid == virus_id]
    }
    if (!is.null(il6_id) && any(qid == il6_id)) {
      f_pred <- ode_df$F[time_index]
      pred[qid == il6_id] <- f_pred[qid == il6_id]
    }
    if (!is.null(wl_id) && any(qid == wl_id)) {
      s_pred <- ode_df$S[time_index]
      pred[qid == wl_id] <- s_pred[qid == wl_id]
    }

    add <- sigma_add_by_id[qid]
    prop <- sigma_prop_by_id[qid]
    variance <- pmax(add^2 + (prop * pred)^2, objective_data$min_variance)
    residual <- scenario_data$value - pred
    nll_point <- 0.5 * (log(variance) + (residual^2) / variance)
    weighted_nll_point <- scenario_data$weight * nll_point

    tibble(
      Quantity = factor(quantity_levels[qid], levels = quantity_levels),
      Scenario = factor(scenario_key, levels = scenario_levels),
      weighted_nll = weighted_nll_point
    ) %>%
      group_by(Quantity, Scenario) %>%
      summarise(weighted_nll = sum(weighted_nll), .groups = "drop")
  })

  bind_rows(block_rows)
}

# -----------------------------------------------------------------------------
# Simulation helper
# -----------------------------------------------------------------------------
solvertype <- "vode"
tols <- 1e-10
tfinal <- 7
dt <- 0.1
output_times <- seq(0, tfinal, by = dt)

run_model_once <- function(params, fixedpars, Y0_vals, simulatorname) {
  params <- params[!is.na(params)]
  fixedpars <- fixedpars[!is.na(fixedpars)]

  params_ode <- params[!grepl("^sigma_", names(params))]
  fixedpars_ode <- fixedpars[!grepl("^sigma_", names(fixedpars))]

  simulate_one <- function(ad0, scenario_label, times) {
    allpars <- c(
      as.list(Y0_vals),
      as.list(params_ode),
      as.list(fixedpars_ode),
      list(
        Ad0 = ad0,
        txstart = 1,
        txinterval = 0.5,
        # Match fitting/simulation schedule: last dose just before day 4
        txend = 3.9,
        tstart = 0,
        tfinal = tfinal,
        dt = dt,
        solvertype = solvertype,
        tols = tols,
        times = times
      )
    )

    odeout <- tryCatch(do.call(simulatorname, allpars), error = identity)
    if (inherits(odeout, "error")) {
      return(odeout)
    }

    ode_df <- as.data.frame(odeout)

    Ct <- ode_df$At / as.numeric(allpars["Vt"])
    fu <- as.numeric(allpars["fmax"]) * Ct / (as.numeric(allpars["f50"]) + Ct)
    Cu <- fu * Ct
    fV <- as.numeric(allpars["Emax_V"]) * Cu / (as.numeric(allpars["C50_V"]) + Cu)
    fF <- as.numeric(allpars["Emax_F"]) * Cu / (as.numeric(allpars["C50_F"]) + Cu)

    vprod <- (1 - fV) * as.numeric(allpars["p"]) * ode_df$I / (1 + as.numeric(allpars["kF"]) * ode_df$F)

    ode_df$Cu <- Cu
    ode_df$fV <- fV
    ode_df$fF <- vprod

    ode_df$Dose <- ad0
    ode_df$Scenario <- factor(scenario_label, levels = scenarios)

    ode_df
  }

  sim_obj_list <- lapply(seq_along(doses), function(i) {
    scenario_key <- as.character(scenarios[i])
    times <- obs_times_by_scenario[[scenario_key]]
    if (is.null(times) || !length(times)) {
      times <- output_times
    }
    simulate_one(doses[i], scenarios[i], times)
  })

  failed <- vapply(sim_obj_list, inherits, logical(1), "error")
  if (any(failed)) {
    first_error <- sim_obj_list[[which(failed)[1]]]
    return(list(error = first_error))
  }

  sim_obj_df <- bind_rows(sim_obj_list)

  objective <- fit_model_function(
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
    simulatorname = simulatorname,
    logfit = 0,
    obs_times_by_scenario = obs_times_by_scenario,
    objective_data = objective_data,
    weight_mode = weight_mode
  )

  sigma_pool <- collect_sigma_pool(params, fixedpars)

  breakdown_nll_fast <- compute_fast_weighted_nll_breakdown(
    sim_obj_df = sim_obj_df,
    sigma_pool = sigma_pool,
    objective_data = objective_data
  )

  pred_long <- build_prediction_long(sim_obj_df, scenario_col = "Scenario", time_col = "time")
  components <- compute_objective_components(
    fitdata,
    pred_long,
    sigma_pool,
    weight_mode = weight_mode
  )
  if (!is.finite(objective)) {
    objective <- simpleError("Objective evaluation failed (non-finite result).")
  }

  if (is.null(breakdown_nll_fast) || !nrow(breakdown_nll_fast)) {
    breakdown_nll_fast <- components$breakdown_nll
  }

  breakdown <- list(
    nll = breakdown_nll_fast,
    ssr = components$breakdown_ssr
  )

  sim_plot_list <- lapply(seq_along(doses), function(i) {
    simulate_one(doses[i], scenarios[i], output_times)
  })

  failed_plot <- vapply(sim_plot_list, inherits, logical(1), "error")
  if (any(failed_plot)) {
    first_error <- sim_plot_list[[which(failed_plot)[1]]]
    return(list(error = first_error))
  }

  sim_plot_df <- bind_rows(sim_plot_list)

  list(sim = sim_plot_df, objective = objective, breakdown = breakdown)
}

# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------
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
        choices = c("Model 1" = "model1", "Model 2" = "model2"),
        selected = "model1"
      ),
      helpText(
        "Adjust the parameters to immediately update the simulation ",
        "and objective function evaluation."
      ),
      br(),
      tabsetPanel(
        id = "param-tabs",
        tabPanel("Fitted parameters", uiOutput("fit_param_inputs")),
        tabPanel(
          "Fixed parameters",
          uiOutput("fixed_param_inputs"),
          uiOutput("fixed_sigma_inputs")
        )
      )
    ),
    mainPanel(
      width = 8,
      h4("Objective function"),
      textOutput("objective_value"),
      br(),
      h4("Objective contributions (by variable and dose)"),
      tableOutput("objective_breakdown"),
      br(),
      h4("SSR contributions (by variable and dose)"),
      tableOutput("objective_ssr"),
      br(),
      plotOutput("fit_plot", height = "900px")
    )
  )
)

# -----------------------------------------------------------------------------
# Server
# -----------------------------------------------------------------------------
server <- function(input, output, session) {
  model_config <- reactive({
    req(input$model_choice)
    model_configs[[input$model_choice]]
  })

  output$fit_param_inputs <- renderUI({
    cfg <- model_config()

    div(
      class = "parameter-grid",
      lapply(cfg$fit_param_names, function(nm) {
        label <- format_input_label(nm, fit_param_labels)
        numericInput(
          inputId = param_input_id(nm),
          label = label,
          value = cfg$fit_defaults[[nm]],
          min = 0,
          step = NA
        )
      })
    )
  })

  output$fixed_param_inputs <- renderUI({
    cfg <- model_config()

    div(
      class = "parameter-grid",
      lapply(cfg$fixed_param_names, function(nm) {
        label <- format_input_label(nm, cfg$fixed_param_labels)
        numericInput(
          inputId = fixed_input_id(nm),
          label = label,
          value = cfg$fixed_defaults[[nm]],
          min = 0,
          step = NA
        )
      })
    )
  })

  output$fixed_sigma_inputs <- renderUI({
    cfg <- model_config()
    sigma_fixed_names <- names(sigma_settings$sigma_fixed)

    div(
      class = "parameter-grid",
      tags$div(class = "control-header", "Measurement noise (fixed)"),
      lapply(sigma_fixed_names, function(nm) {
        label <- format_input_label(nm, sigma_fixed_labels)
        numericInput(
          inputId = fixed_input_id(nm),
          label = label,
          value = cfg$fixed_defaults[[nm]],
          min = 0,
          step = NA
        )
      })
    )
  })

  fitted_values <- reactive({
    cfg <- model_config()
    defaults <- cfg$fit_defaults

    vapply(
      cfg$fit_param_names,
      function(nm) input[[param_input_id(nm)]] %||% defaults[[nm]],
      numeric(1),
      USE.NAMES = TRUE
    )
  })

  fixed_values <- reactive({
    cfg <- model_config()
    defaults <- cfg$fixed_defaults
    fixed_names <- c(cfg$fixed_param_names, names(sigma_settings$sigma_fixed))

    vapply(
      fixed_names,
      function(nm) input[[fixed_input_id(nm)]] %||% defaults[[nm]],
      numeric(1),
      USE.NAMES = TRUE
    )
  })

  model_results <- reactive({
    cfg <- model_config()
    params <- fitted_values()
    fixedpars <- fixed_values()

    run_model_once(params, fixedpars, cfg$Y0, cfg$simulatorname)
  })

  output$objective_value <- renderText({
    res <- model_results()

    if (!is.null(res$error)) {
      return(paste("Simulation failed:", conditionMessage(res$error)))
    }

    if (inherits(res$objective, "error")) {
      return(paste("Objective evaluation failed:", conditionMessage(res$objective)))
    }

    paste0("Current objective value: ", format(res$objective, digits = 6, scientific = TRUE))
  })

  output$objective_breakdown <- renderTable({
    res <- model_results()
    if (is.null(res$breakdown) || is.null(res$breakdown$nll)) {
      return(NULL)
    }

    format_breakdown_table(res$breakdown$nll, "weighted_nll")
  }, digits = 4, rownames = TRUE)

  output$objective_ssr <- renderTable({
    res <- model_results()
    if (is.null(res$breakdown) || is.null(res$breakdown$ssr)) {
      return(NULL)
    }

    format_breakdown_table(res$breakdown$ssr, "weighted_ssr")
  }, digits = 4, rownames = TRUE)

  output$fit_plot <- renderPlot({
    res <- model_results()
    validate(need(is.null(res$error), if (is.null(res$error)) "" else conditionMessage(res$error)))
    validate(need(!inherits(res$objective, "error"),
                  if (inherits(res$objective, "error")) conditionMessage(res$objective) else ""))

    sim_df <- res$sim %>%
      mutate(Dose = factor(Dose, levels = doses))

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

shinyApp(ui, server)
