# -----------------------------------------------------------------------------
# diagnostic-plot-function.R
# -----------------------------------------------------------------------------
# PURPOSE
#   Build residual and predicted-vs-observed diagnostic plots for fitted models.
#   These plots are consistent with the NLL weighting used in the objective.
# -----------------------------------------------------------------------------

library(ggplot2) # Plotting.
library(dplyr)   # Data manipulation.
library(here)    # Project-root-relative paths.

# Centralized objective computations (shared with fitter + app).
source(here::here("code", "analysis-code", "functions", "objective-components-function.R"))

#' Create diagnostic plots for a single bestfit.
#'
#' @param bestfit A single bestfit object (one list element).
#' @param sim_df Time-series simulations from dose-prediction output (one element).
#' @param output_prefix File prefix to use when saving figures.
#' @param figures_dir Directory to save figures in.
#' @param config Plotting config from get_plotting_config().
#' @param jitter_width_combined Jitter width for combined residual plot.
#' @param jitter_width_panel Jitter width for panel residual plots.
#' @return A list of ggplot objects.
#'
create_diagnostic_plots <- function(
  bestfit,
  sim_df,
  output_prefix,
  figures_dir,
  config,
  jitter_width_combined = 0.08,
  jitter_width_panel = 0.08
) {
  # Guard against missing or empty simulation data.
  if (is.null(sim_df) || !nrow(sim_df)) {
    stop("Simulation data is missing; cannot build diagnostics.")
  }

  # Observed data from the bestfit object.
  dat <- bestfit$fitdata

  # Enforce consistent factor ordering.
  dat$Scenario <- factor(dat$Scenario, levels = config$scenario_levels)
  dat$Quantity <- factor(dat$Quantity, levels = config$quantity_levels)

  # Keep only the baseline schedule and experimental doses (0/10/100).
  dose_levels <- as.numeric(names(config$scenario_map))

  # Convert simulation output into long format aligned to observed data times.
  sim_df_use <- sim_df %>%
    filter(Schedule == "s1", Dose %in% dose_levels) %>%
    mutate(
      Scenario = factor(
        config$scenario_map[as.character(Dose)],
        levels = config$scenario_levels
      )
    )

  pred_long <- build_prediction_long(sim_df_use, scenario_col = "Scenario", time_col = "time")

  sigma_pool <- collect_sigma_pool(bestfit$fitpars, bestfit$fixedpars)
  components <- compute_objective_components(dat, pred_long, sigma_pool)
  resid_df <- components$residuals

  if (is.null(resid_df) || !nrow(resid_df)) {
    stop("No residuals were computed; check simulation/data alignment.")
  }

  resid_df <- resid_df %>%
    mutate(
      Quantity = factor(Quantity, levels = config$quantity_levels),
      Scenario = factor(Scenario, levels = config$scenario_levels)
    )

  scen_levels <- levels(resid_df$Scenario)
  shape_vals <- setNames(c(16, 17, 15)[seq_along(scen_levels)], scen_levels)

  # Combined residual plot (all quantities together).
  y_max_comb <- max(abs(resid_df$NLLWeightedResid), na.rm = TRUE)
  if (!is.finite(y_max_comb) || y_max_comb == 0) y_max_comb <- 1
  y_pad_comb <- 0.1 * y_max_comb
  y_lim_comb <- c(-(y_max_comb + y_pad_comb), y_max_comb + y_pad_comb)

  p_combined <- ggplot(
    resid_df,
    aes(x = Day, y = NLLWeightedResid, colour = Quantity, shape = Scenario)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(
      position = position_jitter(width = jitter_width_combined, height = 0),
      alpha = 0.85,
      size = 2
    ) +
    scale_color_manual(
      values = config$var_colors,
      name = "Variable",
      labels = config$var_labs,
      guide = guide_legend(order = 1)
    ) +
    scale_shape_manual(
      values = shape_vals,
      name = "Scenario",
      labels = config$scen_labs,
      guide = guide_legend(order = 2)
    ) +
    scale_y_continuous(limits = y_lim_comb, expand = expansion(mult = 0)) +
    xlab("Time (days)") +
    ylab("NLL-weighted standardized residuals") +
    theme_bw(base_size = 14) +
    theme(legend.position = "top", legend.box = "vertical")

  # Faceted residual plot (one panel per quantity).
  pad_df <- resid_df %>%
    group_by(Quantity) %>%
    summarise(
      ypad = max(abs(NLLWeightedResid), na.rm = TRUE) * 1.1,
      x_min = min(Day, na.rm = TRUE),
      x_max = max(Day, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    transmute(Quantity, ymin = -ypad, ymax = ypad, x_min, x_max)

  p_panel <- ggplot(
    resid_df,
    aes(x = Day, y = NLLWeightedResid, colour = Scenario, shape = Scenario)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(
      position = position_jitter(width = jitter_width_panel, height = 0),
      alpha = 0.85,
      size = 2
    ) +
    geom_blank(data = pad_df, aes(x = x_min, y = ymin), inherit.aes = FALSE) +
    geom_blank(data = pad_df, aes(x = x_max, y = ymax), inherit.aes = FALSE) +
    facet_wrap(
      ~Quantity,
      nrow = 1,
      labeller = labeller(Quantity = config$var_labs),
      scales = "free_y"
    ) +
    scale_color_manual(
      values = config$var_colors,
      name = "Scenario",
      labels = config$scen_labs
    ) +
    scale_shape_manual(values = shape_vals, name = "Scenario", labels = config$scen_labs) +
    xlab("Time (days)") +
    ylab("NLL-weighted standardized residuals") +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "top",
      strip.background = element_blank(),
      strip.text = element_text(size = 12)
    )

  # Predicted vs observed, faceted by quantity.
  gof_df <- resid_df %>%
    transmute(Scenario, Quantity, Observed = Value, Predicted = Predicted)

  lims_df <- gof_df %>%
    group_by(Quantity) %>%
    summarise(
      low = min(c(Observed, Predicted), na.rm = TRUE),
      high = max(c(Observed, Predicted), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      pad = 0.05 * (high - low + ifelse(high - low == 0, 1, 0)),
      low = low - pad,
      high = high + pad
    )

  lims_pts <- bind_rows(
    lims_df %>% transmute(Quantity, x = low, y = low),
    lims_df %>% transmute(Quantity, x = high, y = high)
  )

  p_gof <- ggplot(
    gof_df,
    aes(x = Observed, y = Predicted, colour = Scenario, shape = Scenario)
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.85, size = 2) +
    geom_blank(data = lims_pts, aes(x = x, y = y), inherit.aes = FALSE) +
    facet_wrap(
      ~Quantity,
      nrow = 1,
      labeller = labeller(Quantity = config$var_labs),
      scales = "free"
    ) +
    theme(aspect.ratio = 1) +
    scale_color_manual(values = config$var_colors, name = "Scenario:", labels = config$scen_labs) +
    scale_shape_manual(values = shape_vals, name = "Scenario:", labels = config$scen_labs) +
    xlab("Observed") +
    ylab("Predicted") +
    theme_bw(base_size = 14) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 12),
      legend.position = "top"
    )

  # Save figures to disk.
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  ggsave(
    filename = file.path(figures_dir, paste0(output_prefix, "-residuals-combined.png")),
    plot = p_combined,
    width = 8,
    height = 5,
    units = "in",
    dpi = 300
  )

  ggsave(
    filename = file.path(figures_dir, paste0(output_prefix, "-residuals-panels.png")),
    plot = p_panel,
    width = 8.5,
    height = 5,
    units = "in",
    dpi = 300
  )

  ggsave(
    filename = file.path(figures_dir, paste0(output_prefix, "-pred-vs-obs.png")),
    plot = p_gof,
    width = 8,
    height = 5,
    units = "in",
    dpi = 300
  )

  list(
    residuals_combined = p_combined,
    residuals_panel = p_panel,
    predicted_vs_observed = p_gof
  )
}
