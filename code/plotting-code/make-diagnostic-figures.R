##########################################
# make residual plot for best fit
##########################################
############################################
##  Packages
############################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork) # for combining ggplot objects
library(here) # for file paths
library(grid) # for a custom x-axis label


############################################
##  1.  Read in the fit results
############################################

bestfit_list = readRDS(here::here('results', 'output', 'bestfit.Rds'))


nsamp = 1
# uncomment the line below to generate best fit tables for all samples of the fixed parameters
#nsamp = length(bestfit_list)
for (i in 1:nsamp) {
  bestfit = bestfit_list[[i]]
  dat <- bestfit$fitdata
  sim <- bestfit$simresult

  times_to_keep <- unique(dat$xvals)

  sim_long <- sim %>%
    filter(time %in% times_to_keep) %>% # keep only matching times
    select(time, Scenario, V, F, S) %>% # keep relevant columns
    mutate(
      V = log10(pmax(1, V)) #take log, also adjust for censoring
    ) %>%
    rename(Day = time, LogVirusLoad = V, IL6 = F, Weight = S) %>%
    pivot_longer(
      cols = c(LogVirusLoad, IL6, Weight),
      names_to = "Quantity",
      values_to = "Predicted"
    )

  ###################################
  # 3.  One-to-one join & residuals
  ###################################
  # Maximum value for each quantity (for weighted residuals)
  max_lookup <- dat %>%
    group_by(Quantity) %>%
    summarise(
      ScalingMax = max(Value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      ScalingMax = ifelse(is.finite(ScalingMax), ScalingMax, NA_real_)
    )

  # Sample size by (Scenario, Quantity) pair to match the objective function
  # scaling from `fit_model_function()`.
  n_lookup <- dat %>%
    group_by(Scenario, Quantity) %>%
    summarise(
      SampleSize = dplyr::n(),
      .groups = "drop"
    )

  resid_df <- dat %>% # keep every data row (incl. replicates)
    inner_join(
      sim_long, # exact match on variable, time, scenario
      by = c("Scenario", "Day", "Quantity")
    ) %>%
    mutate(
      Residual = Value - Predicted,
      # facet order
      Quantity = factor(Quantity, levels = c("LogVirusLoad", "IL6", "Weight")),
      Scenario = factor(
        Scenario,
        levels = c("NoTreatment", "PanCytoVir10mg", "PanCytoVir100mg")
      )
    )

  resid_df <- resid_df %>%
    left_join(max_lookup, by = "Quantity") %>%
    left_join(n_lookup, by = c("Scenario", "Quantity")) %>%
    mutate(
      WeightedResidual = dplyr::if_else(
        !is.na(ScalingMax) &
          ScalingMax > 0 &
          !is.na(SampleSize) &
          SampleSize > 0,
        Residual / ScalingMax / sqrt(SampleSize),
        NA_real_
      )
    )

  ###################################
  # 4.  Nice facet labels
  ###################################
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

  ###################################
  # 5.  Build the 1×3 residual plot (all scenarios overlaid in each panel)
  ###################################

  # colour/shape maps (Okabe–Ito palette)
  scen_levels <- levels(resid_df$Scenario)
  col_vals <- setNames(
    c("#0072B2", "#009E73", "#D55E00")[seq_along(scen_levels)],
    scen_levels
  )
  shape_vals <- setNames(c(16, 17, 15)[seq_along(scen_levels)], scen_levels)

  # Helper that builds the residual plot for a selected residual definition.
  build_combined_residual_plot <- function(df, y_label) {
    df <- df %>%
      filter(!is.na(Residual))

    if (nrow(df) == 0) {
      stop(
        "No residual values available for plotting. Check the weighting inputs."
      )
    }

    pad_df <- df %>%
      dplyr::group_by(Quantity) %>%
      dplyr::summarise(
        max_abs = max(abs(Residual), na.rm = TRUE),
        max_abs = ifelse(is.finite(max_abs), max_abs, 0),
        ypad = max_abs + 0.1 * max_abs,
        x_min = min(Day, na.rm = TRUE),
        x_max = max(Day, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::transmute(Quantity, ymin = -ypad, ymax = ypad, x_min, x_max)

    base_plot <- ggplot(
      df,
      aes(x = Day, y = Residual, colour = Scenario, shape = Scenario)
    ) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_point(alpha = 0.75, size = 2) +
      # --- IMPORTANT: don't inherit colour/shape (no Scenario column in pad_df)
      geom_blank(data = pad_df, aes(x = x_min, y = ymin), inherit.aes = FALSE) +
      geom_blank(data = pad_df, aes(x = x_max, y = ymax), inherit.aes = FALSE) +
      facet_wrap(
        ~Quantity,
        nrow = 1,
        labeller = labeller(Quantity = var_labs),
        scales = "free_y"
      ) +
      scale_y_continuous(expand = expansion(mult = 0)) +
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
      ylab(y_label) +
      theme_bw(base_size = 14) +
      theme(
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "top"
      )

    blank <- plot_spacer()

    full_plot <-
      base_plot /
      blank + # stack blank row below
      plot_layout(heights = c(1, 0.03)) & # tiny height for spacer
      theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5))

    return(full_plot)
  }

  ###################################
  # 6.  Unweighted residual plot (original definition)
  ###################################
  unweighted_plot <- build_combined_residual_plot(resid_df, "Residual")

  print(unweighted_plot)
  grid.text("Time (days)", y = unit(0.02, "npc"), gp = gpar(fontsize = 14))

  figname = paste0('residuals-combined', i, '.png')
  png(
    here::here('results', 'figures', figname),
    width = 8,
    height = 5,
    units = 'in',
    res = 300
  )
  print(unweighted_plot)
  dev.off()

  ###################################
  # 7.  Weighted residual plot (objective-function scale)
  ###################################
  weighted_resid_df <- resid_df %>%
    mutate(Residual = WeightedResidual)

  weighted_plot <- build_combined_residual_plot(
    weighted_resid_df,
    "Weighted Residual"
  )

  print(weighted_plot)
  grid.text("Time (days)", y = unit(0.02, "npc"), gp = gpar(fontsize = 14))

  figname_weighted = paste0('residuals-weighted-combined', i, '.png')
  png(
    here::here('results', 'figures', figname_weighted),
    width = 8,
    height = 5,
    units = 'in',
    res = 300
  )
  print(weighted_plot)
  dev.off()

  ###################################
  # 7b. Weighted residual plot — SINGLE PANEL
  #     color = Quantity (virus, IL-6, weight)
  #     shape = Scenario (dose level)
  ###################################

  # Okabe–Ito colors mapped to variables
  var_colors <- c(
    LogVirusLoad = "#0072B2", # blue
    IL6 = "#D55E00", # vermillion
    Weight = "#009E73" # bluish-green
  )

  # Ensure factors are ordered nicely
  weighted_single_df <- weighted_resid_df %>%
    filter(is.finite(Residual)) %>%
    mutate(
      Quantity = factor(Quantity, levels = c("LogVirusLoad", "IL6", "Weight")),
      Scenario = factor(Scenario, levels = levels(resid_df$Scenario))
    )

  # Symmetric y-limits around zero (with a little padding)
  y_max <- max(abs(weighted_single_df$Residual), na.rm = TRUE)
  if (!is.finite(y_max) || y_max == 0) {
    y_max <- 1
  }
  y_pad <- 0.1 * y_max
  y_limits <- c(-(y_max + y_pad), y_max + y_pad)

  weighted_singlepanel_plot <- ggplot(
    weighted_single_df,
    aes(x = Day, y = Residual, colour = Quantity, shape = Scenario)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.85, size = 2) +
    scale_color_manual(
      values = var_colors,
      name = "Variable:",
      labels = var_labs
    ) +
    scale_shape_manual(
      values = shape_vals, # uses the shapes you defined earlier for scenarios
      name = "Dose level:",
      labels = scen_labs
    ) +
    scale_y_continuous(limits = y_limits, expand = expansion(mult = 0)) +
    xlab("Time (days)") +
    ylab("Weighted Residual") +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "top",
      strip.background = element_blank()
    )

  print(weighted_singlepanel_plot)

  figname_weighted_single <- paste0('residuals-weighted-singlepanel', i, '.png')
  png(
    here::here('results', 'figures', figname_weighted_single),
    width = 7.5,
    height = 5,
    units = 'in',
    res = 300
  )
  print(weighted_singlepanel_plot)
  dev.off()

  ###################################
  # 7. Predicted vs Observed GOF plot
  ###################################
  # Reuse the joined residual dataframe to get Observed & Predicted
  gof_df <- resid_df %>%
    transmute(
      Scenario,
      Quantity,
      Observed = Value,
      Predicted = Predicted
    )

  # Per-facet limits so x and y share identical ranges (with small padding)
  lims_df <- gof_df %>%
    dplyr::group_by(Quantity) %>%
    dplyr::summarise(
      low = min(c(Observed, Predicted), na.rm = TRUE),
      high = max(c(Observed, Predicted), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      pad = 0.05 * (high - low + ifelse(high - low == 0, 1, 0)), # guard zero-range
      low = low - pad,
      high = high + pad
    )

  # Invisible points to enforce identical x/y ranges per facet
  lims_pts <- dplyr::bind_rows(
    lims_df %>% dplyr::transmute(Quantity, x = low, y = low),
    lims_df %>% dplyr::transmute(Quantity, x = high, y = high)
  )

  gof_p <- ggplot(
    gof_df,
    aes(x = Observed, y = Predicted, colour = Scenario, shape = Scenario)
  ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      colour = "grey50"
    ) +
    geom_point(alpha = 0.75, size = 2) +
    geom_blank(data = lims_pts, aes(x = x, y = y), inherit.aes = FALSE) +
    facet_wrap(
      ~Quantity,
      nrow = 1,
      labeller = labeller(Quantity = var_labs),
      scales = "free"
    ) +
    theme(aspect.ratio = 1) +
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
    xlab("Observed") +
    ylab("Predicted") +
    theme_bw(base_size = 14) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 12),
      legend.position = "top"
    )

  # Save GOF figure
  figname_gof <- paste0('pred-vs-obs-combined', i, '.png')
  png(
    here::here('results', 'figures', figname_gof),
    width = 8,
    height = 5,
    units = 'in',
    res = 300
  )
  print(gof_p)
  dev.off()
} #end loop over all samples
