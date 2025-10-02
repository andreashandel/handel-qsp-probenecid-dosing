###############################################################################
# make-npde-figures.R
#
# This plotting script generates two diagnostic figures that visualise the
# Normalized Prediction Distribution Errors (NPDE) computed by the
# `run-npde-simulations.R` workflow:
#   * Faceted histograms that compare the empirical NPDE distribution for each
#     outcome variable to the standard normal density.
#   * Quantile-Quantile (QQ) plots that contrast the empirical quantiles with
#     the theoretical quantiles of a standard normal distribution.
#
# The figures are saved in `results/figures/` and are intended for inclusion in
# reports or presentations.  The code is extensively documented so that readers
# with minimal R or statistical background can understand each step.
###############################################################################

############################################
## 1. Load libraries and NPDE results
############################################

library(here)      # platform-independent file paths
library(tidyverse) # includes ggplot2, dplyr, readr, etc.

results_path <- here::here("results", "output", "npde-results.Rds")
if (!file.exists(results_path)) {
  stop("The file 'results/output/npde-results.Rds' was not found. Run 'run-npde-simulations.R' first.")
}

npde_results <- readRDS(results_path)
npde_table <- npde_results$npde_table

# For clarity in the plots we attach human-readable labels to each quantity.
quantity_labels <- c(
  LogVirusLoad = "Log Virus Load",
  IL6 = "IL-6",
  Weight = "Weight"
)

npde_table <- npde_table %>%
  dplyr::mutate(
    Quantity = factor(Quantity, levels = names(quantity_labels), labels = quantity_labels)
  )

############################################
## 2. Histogram diagnostics
############################################

histogram_plot <- ggplot(npde_table, aes(x = NPDE)) +
  geom_histogram(
    bins = 15,
    colour = "white",
    fill = "#4C72B0",
    alpha = 0.75
  ) +
  stat_function(
    fun = dnorm,
    args = list(mean = 0, sd = 1),
    colour = "#DD8452",
    linewidth = 1
  ) +
  facet_wrap(~Quantity, scales = "free_y") +
  labs(
    title = "NPDE Histograms",
    subtitle = "The orange curve shows the standard normal density expected under a well-calibrated model.",
    x = "NPDE",
    y = "Count"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11)
  )

############################################
## 3. QQ-plot diagnostics
############################################

qq_plot <- ggplot(npde_table, aes(sample = NPDE)) +
  stat_qq(colour = "#4C72B0", alpha = 0.8) +
  stat_qq_line(colour = "#DD8452") +
  facet_wrap(~Quantity) +
  labs(
    title = "NPDE Quantile-Quantile Plots",
    subtitle = "Points should follow the orange reference line if NPDE values are normally distributed.",
    x = "Theoretical Quantiles",
    y = "Empirical Quantiles"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11)
  )

############################################
## 4. Save figures to disk
############################################

fig_dir <- here::here("results", "figures")
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

hist_path <- file.path(fig_dir, "npde-histograms.png")
qq_path <- file.path(fig_dir, "npde-qqplots.png")

ggplot2::ggsave(hist_path, plot = histogram_plot, width = 8, height = 4.5, dpi = 300)
ggplot2::ggsave(qq_path, plot = qq_plot, width = 8, height = 4.5, dpi = 300)

message("Saved NPDE histogram to ", hist_path)
message("Saved NPDE QQ plot to ", qq_path)

###############################################################################
# End of script
###############################################################################
