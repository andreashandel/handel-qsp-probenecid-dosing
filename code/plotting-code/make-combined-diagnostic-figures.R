##########################################
# make residual plot for best fit
##########################################
############################################
##  Packages
############################################
library(tidyverse) # dplyr, tidyr, ggplot2, readr, …
library(patchwork) # for combining ggplot objects
library(here) # for file paths
library(grid) # for a custom x-axis label


############################################
##  1.  Read in the fit results
############################################

bestfit_list = readRDS(here::here('results', 'output', 'bestfit.Rds'))


nsamp = 1
# uncomment the line below to generate best fit tables for all samples of the fixed parameters
nsamp = length(bestfit_list)
for (i in 1:nsamp) {
  bestfit = bestfit_list[[i]]
  dat <- bestfit$fitdata
  sim <- bestfit$simresult

  times_to_keep <- unique(dat$xvals)

  sim_long <- sim %>%
    filter(time %in% times_to_keep) %>% # keep only matching times
    select(time, Scenario, V, F, S) %>% # keep relevant columns
    mutate(
      V = log10(V + 1e-12)
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

  # symmetric padding data (as before)
  pad_df <- resid_df %>%
    dplyr::group_by(Quantity) %>%
    dplyr::summarise(
      ypad = max(abs(Residual), na.rm = TRUE) +
        0.1 * max(abs(Residual), na.rm = TRUE),
      x_min = min(Day, na.rm = TRUE),
      x_max = max(Day, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::transmute(Quantity, ymin = -ypad, ymax = ypad, x_min, x_max)

  p <- ggplot(
    resid_df,
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
    ylab("Residual") +
    theme_bw(base_size = 14) +
    theme(
      axis.title.x = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 12),
      legend.position = "top"
    )

  ###################################
  # 6.  Add a single centred x-axis
  ###################################
  #   Patchwork trick: stack a blank
  #   spacer row, then draw the label
  ###################################
  blank <- plot_spacer()

  full_plot <-
    p /
    blank + # stack blank row below
    plot_layout(heights = c(1, 0.03)) & # tiny height for spacer
    theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5))

  # draw the plot and then the label
  print(full_plot)
  grid.text("Time (days)", y = unit(0.02, "npc"), gp = gpar(fontsize = 14))

  # save to png file
  figname = paste0('residuals-combined', i, '.png')
  png(
    here::here('results', 'figures', figname),
    width = 8,
    height = 5,
    units = 'in',
    res = 300
  )
  print(full_plot)
  dev.off()
}
