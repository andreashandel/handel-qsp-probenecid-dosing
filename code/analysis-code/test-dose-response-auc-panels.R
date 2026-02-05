# -----------------------------------------------------------------------------
# test-dose-response-auc-panels.R
# -----------------------------------------------------------------------------
# PURPOSE
#   Testing utility: for each fixed-parameter sample, plot dose vs total virus
#   load (AUC) for the baseline treatment schedule only. Stand-alone script
#   that reads the existing dose-response results file.
# -----------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(here)

# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------
model_choice <- "model1" # "model1" or "model2"

# Input file from run-dose-predictions.R
results_file <- here::here(
  "results",
  "output",
  paste0(model_choice, "-dose-response-results.Rds")
)

# Output figure
output_file <- here::here(
  "results",
  "figures",
  paste0(model_choice, "-dose-response-auc-panels-test.png")
)

# -----------------------------------------------------------------------------
# Time-series settings (virus load only, baseline schedule)
# -----------------------------------------------------------------------------
sample_index <- 6
ts_doses <- c(0, 0.1, 0.3, 1, 3, 10, 30, 900, 950, 1000, 1050, 1100)

solvertype <- "vode"
tols <- 1e-9
dt <- 0.01
tfinal <- 7

ts_output_file <- here::here(
  "results",
  "figures",
  paste0(model_choice, "-timeseries-virusload-panels-test-sample", sample_index, ".png")
)

# -----------------------------------------------------------------------------
# Time-series settings (all variables for one sample + one dose)
# -----------------------------------------------------------------------------
single_sample_index <- 6
single_dose <- 1100

single_ts_output_file <- here::here(
  "results",
  "figures",
  paste0(
    model_choice,
    "-timeseries-allvars-test-sample",
    single_sample_index,
    "-dose",
    format(single_dose, trim = TRUE, scientific = FALSE),
    ".png"
  )
)

# -----------------------------------------------------------------------------
# Load results
# -----------------------------------------------------------------------------
if (!file.exists(results_file)) {
  stop("Dose-response results not found: ", results_file)
}

results_list <- readRDS(results_file)

# -----------------------------------------------------------------------------
# Build plotting data (baseline schedule only)
# -----------------------------------------------------------------------------
plot_df <- bind_rows(lapply(seq_along(results_list), function(i) {
  res <- results_list[[i]]
  df <- res$all_results_df
  if (is.null(df) || !nrow(df)) {
    return(NULL)
  }
  df %>%
    filter(Schedule == "s1") %>%
    mutate(Sample = paste0("sample_", i)) %>%
    select(Sample, Dose, AUCV)
}))

if (is.null(plot_df) || !nrow(plot_df)) {
  stop("No baseline AUC data found in results_list.")
}

plot_df <- plot_df %>%
  mutate(Sample = factor(Sample, levels = unique(Sample))) %>%
  arrange(Sample, Dose)

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
p <- ggplot(plot_df, aes(x = Dose, y = AUCV)) +
  geom_line(linewidth = 0.6, colour = "#0072B2") +
  scale_x_log10() +
  labs(
    x = "Dose (mg/kg)",
    y = "Total viral load (AUC)",
    title = "Baseline dose-response by fixed-parameter sample"
  ) +
  facet_wrap(~Sample, scales = "free_y") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 9),
    plot.title = element_text(size = 14)
  )

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
ggsave(output_file, p, width = 12, height = 8, units = "in", dpi = 300)

message("Saved test dose-response panel plot to: ", output_file)

# -----------------------------------------------------------------------------
# Time-series plots: virus load panels for selected doses (baseline schedule)
# -----------------------------------------------------------------------------
bestfit_file <- here::here(
  "results",
  "output",
  paste0(model_choice, "-bestfit-sample.Rds")
)

if (!file.exists(bestfit_file)) {
  stop("Bestfit file not found: ", bestfit_file)
}

# Load simulator (model1 or model2)
source(here::here("code", "analysis-code", "functions", "model1-simulator-function.R"))
source(here::here("code", "analysis-code", "functions", "model2-simulator-function.R"))
simulatorname <- if (model_choice == "model1") model1_simulator else model2_simulator
source(here::here("code", "plotting-code", "functions", "timeseries-plot-function.R"))

bestfit_list <- readRDS(bestfit_file)
if (sample_index > length(bestfit_list)) {
  stop("sample_index exceeds available bestfit samples.")
}

bestfit <- bestfit_list[[sample_index]]

fit_sigmas <- grepl("^sigma_(add|prop)_", names(bestfit$fitpars))
fitpars_ode <- bestfit$fitpars[!fit_sigmas]

fixed_sigmas <- grepl("^sigma_(add|prop)_", names(bestfit$fixedpars))
fixedpars_ode <- bestfit$fixedpars[!fixed_sigmas]

Y0 <- bestfit$Y0

simulate_one <- function(ad0) {
  allpars <- c(
    as.list(Y0),
    as.list(fitpars_ode),
    as.list(fixedpars_ode),
    list(
      Ad0 = ad0,
      txstart = 1,
      txinterval = 0.5,
      txend = 3.9,
      tstart = 0,
      tfinal = tfinal,
      dt = dt,
      solvertype = solvertype,
      tols = tols
    )
  )

  odeout <- tryCatch(
    do.call(simulatorname, allpars),
    error = function(e) e
  )
  if (inherits(odeout, "error")) {
    warning("Simulation failed for dose ", ad0, ": ", conditionMessage(odeout))
    return(NULL)
  }

  ode_df <- as.data.frame(odeout)
  ode_df$Dose <- ad0
  ode_df
}

ts_list <- lapply(ts_doses, simulate_one)
ts_df <- bind_rows(ts_list)

if (!nrow(ts_df)) {
  stop("No time-series data produced for requested doses.")
}

ts_plot_df <- bind_rows(
  ts_df %>% mutate(DoseFacet = as.character(Dose)),
  ts_df %>% mutate(DoseFacet = "All")
)
ts_plot_df$DoseFacet <- factor(
  ts_plot_df$DoseFacet,
  levels = c("All", as.character(ts_doses))
)

p_ts <- ggplot(ts_plot_df, aes(x = time, y = V, group = Dose)) +
  geom_line(linewidth = 0.6, colour = "#0072B2") +
  facet_wrap(~DoseFacet, scales = "fixed") +
  scale_y_log10() +
  labs(
    x = "Time (days)",
    y = "Virus load (log scale)",
    title = paste0("Baseline virus load trajectories (sample ", sample_index, ")")
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 9),
    plot.title = element_text(size = 14)
  )

ggsave(ts_output_file, p_ts, width = 12, height = 8, units = "in", dpi = 300)
message("Saved test virus-load time-series panels to: ", ts_output_file)

# -----------------------------------------------------------------------------
# Time-series plot: all variables for one sample + one dose (baseline schedule)
# -----------------------------------------------------------------------------
if (single_sample_index > length(bestfit_list)) {
  stop("single_sample_index exceeds available bestfit samples.")
}

bestfit_single <- bestfit_list[[single_sample_index]]

simulate_single_dose <- function(bestfit_obj, ad0) {
  fit_sigmas <- grepl("^sigma_(add|prop)_", names(bestfit_obj$fitpars))
  fitpars_ode <- bestfit_obj$fitpars[!fit_sigmas]

  fixed_sigmas <- grepl("^sigma_(add|prop)_", names(bestfit_obj$fixedpars))
  fixedpars_ode <- bestfit_obj$fixedpars[!fixed_sigmas]

  allpars <- c(
    as.list(bestfit_obj$Y0),
    as.list(fitpars_ode),
    as.list(fixedpars_ode),
    list(
      Ad0 = ad0,
      txstart = 1,
      txinterval = 0.5,
      txend = 3.9,
      tstart = 0,
      tfinal = tfinal,
      dt = dt,
      solvertype = solvertype,
      tols = tols
    )
  )

  odeout <- tryCatch(
    do.call(simulatorname, allpars),
    error = function(e) e
  )
  if (inherits(odeout, "error")) {
    stop("Simulation failed for dose ", ad0, ": ", conditionMessage(odeout))
  }

  ode_df <- as.data.frame(odeout)
  ode_df$Dose <- ad0
  ode_df$Schedule <- "s1"
  ode_df
}

single_sim <- simulate_single_dose(bestfit_single, single_dose)

single_data <- bestfit_single$fitdata %>%
  filter(Dose == single_dose)

format_dose_label <- function(dose) {
  if (isTRUE(all.equal(dose, 0))) {
    return("no drug")
  }
  paste0(format(dose, trim = TRUE, scientific = FALSE), " mg/kg")
}

single_plot <- plot_timeseries(
  data = single_data,
  modelfit = single_sim,
  tmax = tfinal,
  dose_levels = single_dose,
  dose_levels_labels = format_dose_label(single_dose),
  x_jitter = 0.3
)

ggsave(single_ts_output_file, single_plot, width = 8, height = 5, units = "in", dpi = 300)
message("Saved single-dose all-variable time-series plot to: ", single_ts_output_file)
