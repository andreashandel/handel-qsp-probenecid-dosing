##########################################
# make-vpc-figures.R
##########################################
#
# Purpose --------------------------------------------------------------------
# This script creates visual predictive check (VPC) figures for the calibrated
# quantitative systems pharmacology (QSP) model.  The VPC compares the
# experimental observations with the posterior predictive distribution that was
# generated during the NPDE workflow (`run-npde-simulations.R`).  The plot shows
# how well the model reproduces the magnitude and variability of the data across
# all scenarios (No treatment, 10 mg/kg, 100 mg/kg) and outcome variables (virus
# load, IL-6, and body weight).
#
# Overview --------------------------------------------------------------------
# 1. Load the NPDE results, which contain both the observed data and the
#    posterior predictive draws for every observation time.
# 2. Summarise the predictive draws into point-wise quantiles (5th, 50th, and
#    95th percentiles) for each scenario, quantity, and time point.
# 3. Combine the summaries with the observed measurements to produce a faceted
#    VPC figure that overlays the observations on the predictive intervals.
# 4. Save the final figure to `results/figures/vpc-all.png` so that it is readily
#    available for reports and manuscripts.
#
# The code favours clarity over brevity; explicit intermediate variables and
# descriptive comments are used throughout to document each processing step.
##########################################

############################################
## Packages
############################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(here)

############################################
## 1. Load NPDE results and perform safety checks
############################################

npde_path <- here::here("results", "output", "npde-results.Rds")
if (!file.exists(npde_path)) {
  stop(
    "The file 'results/output/npde-results.Rds' is missing. Run 'run-npde-simulations.R' before generating VPC figures."
  )
}

npde_results <- readRDS(npde_path)

observed_data <- npde_results$observed_data
prediction_draws <- npde_results$prediction_draws

if (is.null(observed_data) || nrow(observed_data) == 0) {
  stop("No observed data found in the NPDE results object.")
}

if (is.null(prediction_draws) || nrow(prediction_draws) == 0) {
  stop("No posterior predictions available for the VPC. Check the NPDE simulations.")
}

# Ensure consistent factor ordering for scenarios and measured quantities.  The
# NPDE helpers already enforce this, but we restate the logic here to make the
# plotting script self-contained and explicit.
scenario_levels <- c("NoTreatment", "PanCytoVir10mg", "PanCytoVir100mg")
quantity_levels <- c("LogVirusLoad", "IL6", "Weight")

observed_data <- observed_data %>%
  mutate(
    Scenario = factor(Scenario, levels = scenario_levels, ordered = TRUE),
    Quantity = factor(Quantity, levels = quantity_levels, ordered = TRUE)
  )

prediction_draws <- prediction_draws %>%
  mutate(
    Scenario = factor(Scenario, levels = scenario_levels, ordered = TRUE),
    Quantity = factor(Quantity, levels = quantity_levels, ordered = TRUE)
  )

############################################
## 2. Summarise posterior predictive distribution
############################################

prediction_summary <- prediction_draws %>%
  group_by(Scenario, Quantity, Day) %>%
  summarise(
    lower = stats::quantile(Prediction, probs = 0.05, na.rm = TRUE),
    median = stats::quantile(Prediction, probs = 0.5, na.rm = TRUE),
    upper = stats::quantile(Prediction, probs = 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(median))

############################################
## 3. Create the VPC figure
############################################

# Friendly facet labels reused in other plotting scripts.
scen_labs <- c(
  NoTreatment = "No Treatment",
  PanCytoVir10mg = "10 mg/kg",
  PanCytoVir100mg = "100 mg/kg"
)
var_labs <- c(
  LogVirusLoad = "Log Virus Load",
  IL6 = "IL-6",
  Weight = "Weight"
)

# Palette and shapes chosen to match existing diagnostic figures (Okabe–Ito).
col_vals <- setNames(
  c("#0072B2", "#009E73", "#D55E00")[seq_along(scenario_levels)],
  scenario_levels
)
shape_vals <- setNames(c(16, 17, 15)[seq_along(scenario_levels)], scenario_levels)

vpc_plot <- ggplot() +
  geom_ribbon(
    data = prediction_summary,
    aes(x = Day, ymin = lower, ymax = upper, fill = Scenario),
    alpha = 0.25,
    colour = NA
  ) +
  geom_line(
    data = prediction_summary,
    aes(x = Day, y = median, colour = Scenario),
    size = 1
  ) +
  geom_point(
    data = observed_data,
    aes(x = Day, y = Value, colour = Scenario, shape = Scenario),
    size = 2,
    alpha = 0.8
  ) +
  facet_grid(
    rows = vars(Quantity),
    cols = vars(Scenario),
    labeller = labeller(Quantity = var_labs, Scenario = scen_labs),
    scales = "free_y"
  ) +
  scale_color_manual(
    values = col_vals,
    name = "Scenario:",
    labels = scen_labs
  ) +
  scale_shape_manual(
    values = shape_vals,
    name = "Scenario:",
    labels = scen_labs
  ) +
  scale_fill_manual(values = col_vals, guide = "none") +
  xlab("Time (days)") +
  ylab("Observed and Predicted Values") +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    legend.position = "top"
  )

############################################
## 4. Save the figure
############################################

fig_dir <- here::here("results", "figures")
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

png(
  filename = file.path(fig_dir, "vpc-all.png"),
  width = 10,
  height = 8,
  units = "in",
  res = 300
)
print(vpc_plot)
dev.off()

##########################################
# End of script
##########################################
