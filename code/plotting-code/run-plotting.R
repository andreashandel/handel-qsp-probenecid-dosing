# -----------------------------------------------------------------------------
# run-plotting.R
# -----------------------------------------------------------------------------
# PURPOSE
#   Single entry point to generate all plots and tables from the plotting-code
#   folder. Toggle outputs on/off with the flags below.
# -----------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)
library(here)
library(gt)

# Plotting helpers
source(here::here("code", "plotting-code", "functions", "timeseries-plot-function.R"))
source(here::here("code", "plotting-code", "functions", "dose-response-plot-function.R"))
source(here::here("code", "plotting-code", "functions", "diagnostic-plot-function.R"))
source(here::here("code", "plotting-code", "functions", "plotting-config-function.R"))

# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------
model_choice <- "model1" # "model1" or "model2"
nsamp <- 1               # How many samples to plot (1 = baseline only)

# Toggle outputs
make_timeseries_figures <- TRUE
make_diagnostic_figures <- TRUE
make_dose_response_figures <- TRUE
make_parameter_tables <- TRUE

# Dose-response sample display:
#   - "band": show uncertainty bands (current default)
#   - "lines": show each sample as a thin line; baseline remains thick
#dose_response_sample_display <- "band"
dose_response_sample_display <- "lines"


# Custom x-axis ranges (dose in mg/kg) for dose-response figures
x_axis_ranges <- list(
  baseline = c(1e-2, 1e5),
  txstart = c(1e-2, 1e5),
  txinterval = c(1e-2, 1e5)
)

# -----------------------------------------------------------------------------
# Execution timing
# -----------------------------------------------------------------------------
start_time <- Sys.time()
message("Starting run-plotting for model: ", model_choice)

# -----------------------------------------------------------------------------
# Load inputs shared across multiple plot types
# -----------------------------------------------------------------------------
results_file <- here::here("results", "output", paste0(model_choice, "-dose-response-results.Rds"))
bestfit_file <- here::here("results", "output", paste0(model_choice, "-bestfit-sample.Rds"))

if (!file.exists(results_file)) {
  stop("Dose-response results not found: ", results_file)
}
if (!file.exists(bestfit_file)) {
  stop("Bestfit file not found: ", bestfit_file)
}

results_list <- readRDS(results_file)
bestfit_list <- readRDS(bestfit_file)

bestfit_used <- attr(results_list, "bestfit_file")
if (!is.null(bestfit_used) && !identical(bestfit_used, bestfit_file)) {
  warning(
    "Dose-response results were generated from a different bestfit file (",
    bestfit_used,
    ") than the one currently loaded for plotting (",
    bestfit_file,
    ")."
  )
}

nsamp <- min(nsamp, length(bestfit_list))

fig_dir <- here::here("results", "figures")
table_dir <- here::here("results", "tables")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

config <- get_plotting_config()

# -----------------------------------------------------------------------------
# Time-series figures
# -----------------------------------------------------------------------------
if (make_timeseries_figures) {
  message("Generating time-series figures...")

  timeseries_doses <- attr(results_list, "ts_doses")

  format_dose_label <- function(dose) {
    if (isTRUE(all.equal(dose, 0))) {
      return("no drug")
    }
    paste0(format(dose, trim = TRUE, scientific = FALSE), " mg/kg")
  }

  supp_dose_labels <- vapply(timeseries_doses, format_dose_label, character(1))
  df_list <- lapply(results_list, `[[`, "timeseries_df")
  all_ts <- df_list[[1]]

  schedule_files <- c(
    s1 = paste0(model_choice, "-timeseries-baseline.png"),
    s2 = paste0(model_choice, "-timeseries-d2tx.png"),
    s3 = paste0(model_choice, "-timeseries-d3tx.png"),
    s4 = paste0(model_choice, "-timeseries-dailytx.png"),
    s5 = paste0(model_choice, "-timeseries-singletx.png")
  )

  for (schedule_id in names(schedule_files)) {
    schedule_df <- all_ts %>%
      filter(Schedule == schedule_id)

    ts_plot <- plot_timeseries(
      data = NULL,
      modelfit = schedule_df,
      tmax = 7,
      dose_levels = timeseries_doses,
      dose_levels_labels = supp_dose_labels
    )

    ggsave(
      here::here(fig_dir, schedule_files[[schedule_id]]),
      ts_plot,
      width = 8,
      height = 5,
      units = "in",
      dpi = 300
    )
  }

  bestfit_doses <- c(0, 10, 100)
  bestfit_dose_labels <- c("no drug", "10 mg/kg", "100 mg/kg")

  for (i in seq_len(nsamp)) {
    bestfit <- bestfit_list[[i]]
    all_ts <- df_list[[i]] %>%
      filter(Schedule == "s1", Dose %in% bestfit_doses)

    bestfit_plots <- plot_timeseries(
      data = bestfit$fitdata,
      modelfit = all_ts,
      tmax = 7,
      dose_levels = bestfit_doses,
      dose_levels_labels = bestfit_dose_labels,
      x_jitter = 0.3
    )

    ggsave(
      here::here(fig_dir, paste0(model_choice, "-bestfit", i, ".png")),
      bestfit_plots,
      width = 8,
      height = 5,
      units = "in",
      dpi = 300
    )
  }
}

# -----------------------------------------------------------------------------
# Diagnostic figures
# -----------------------------------------------------------------------------
if (make_diagnostic_figures) {
  message("Generating diagnostic figures...")

  for (i in seq_len(nsamp)) {
    bestfit <- bestfit_list[[i]]
    sim_df <- results_list[[i]]$timeseries_df

    create_diagnostic_plots(
      bestfit = bestfit,
      sim_df = sim_df,
      output_prefix = paste0(model_choice, "-sample", i),
      figures_dir = fig_dir,
      config = config
    )
  }
}

# -----------------------------------------------------------------------------
# Dose-response figures
# -----------------------------------------------------------------------------
if (make_dose_response_figures) {
  message("Generating dose-response figures...")

  fig1 <- plot_outcomes(
    results_list,
    scenarios = "baseline",
    x_limits = x_axis_ranges$baseline,
    sample_display = dose_response_sample_display
  )

  fig2 <- plot_outcomes(
    results_list,
    scenarios = c("baseline", "d2 start", "d3 start"),
    x_limits = x_axis_ranges$txstart,
    sample_display = dose_response_sample_display
  )

  fig3 <- plot_outcomes(
    results_list,
    scenarios = c("baseline", "daily tx", "single tx"),
    x_limits = x_axis_ranges$txinterval,
    sample_display = dose_response_sample_display
  )

  ggsave(here::here(fig_dir, paste0(model_choice, "-dose-response-baseline.png")), fig1, width = 12, height = 4)
  ggsave(here::here(fig_dir, paste0(model_choice, "-dose-response-txstart.png")), fig2, width = 12, height = 4)
  ggsave(here::here(fig_dir, paste0(model_choice, "-dose-response-txinterval.png")), fig3, width = 12, height = 4)
}

# -----------------------------------------------------------------------------
# Parameter tables (gt) + PNG export
# -----------------------------------------------------------------------------
if (make_parameter_tables) {
  message("Generating parameter tables...")

  for (i in seq_len(nsamp)) {
    bestfit <- bestfit_list[[i]]
    parnames <- bestfit$fitparnames
    parvalues <- bestfit$fitpars
    parlabels <- bestfit$parlabels

    tab <- data.frame(parnames, parvalues, parlabels)

    gttab <- gt(tab) %>%
      cols_label(
        parnames = "Parameter",
        parvalues = "Value",
        parlabels = "Label"
      )

    saveRDS(gttab, file = here::here(table_dir, paste0(model_choice, "-parametertab", i, ".rds")))

    gt::gtsave(
      data = gttab,
      filename = here::here(fig_dir, paste0(model_choice, "-parametertab", i, ".png"))
    )
  }
}

# -----------------------------------------------------------------------------
# Final timing report
# -----------------------------------------------------------------------------
end_time <- Sys.time()
elapsed_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
message("run-plotting finished for model: ", model_choice)
message(sprintf("Total elapsed time: %.2f minutes", elapsed_minutes))
