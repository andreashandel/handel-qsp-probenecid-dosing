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
      V = log10(pmax(1, V))
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
  # Maximum value for each quantity (used in the weighted residual definition)
  max_lookup <- dat %>%
    group_by(Quantity) %>%
    summarise(
      ScalingMax = max(Value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      ScalingMax = ifelse(is.finite(ScalingMax), ScalingMax, NA_real_)
    )

  # Sample size for each (Scenario, Quantity) pair, mirroring the weights used
  # in the least-squares objective function inside `fit_model_function()`.
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
        !is.na(ScalingMax) & ScalingMax > 0 & !is.na(SampleSize) & SampleSize > 0,
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
  # 5.  Helper to build the 3×3 residual plot
  ###################################
  build_residual_plot <- function(df, y_label) {
    df <- df %>%
      filter(!is.na(Residual))

    if (nrow(df) == 0) {
      stop("No residual values available for plotting. Check the weighting inputs.")
    }

    base_plot <- ggplot(df, aes(x = Day, y = Residual)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_point(alpha = 0.65, size = 2) +
      facet_grid(
        rows = vars(Quantity),
        cols = vars(Scenario),
        labeller = labeller(Quantity = var_labs, Scenario = scen_labs),
        scales = "free_y" # tighter scale for IL-6 row
      ) +
      ylab(y_label) +
      theme_bw(base_size = 14) +
      theme(
        axis.title.x = element_blank(), # master label added later
        strip.background = element_blank(),
        strip.text = element_text(size = 12)
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
  # 6.  Unweighted residual plot (original)
  ###################################
  unweighted_plot <- build_residual_plot(resid_df, "Residual")

  print(unweighted_plot)
  grid.text("Time (days)", y = unit(0.02, "npc"), gp = gpar(fontsize = 14))

  figname = paste0('residuals', i, '.png')
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
  # 7.  Weighted residual plot (matches the objective function weighting)
  ###################################
  weighted_resid_df <- resid_df %>%
    mutate(Residual = WeightedResidual)

  weighted_plot <- build_residual_plot(weighted_resid_df, "Weighted Residual")

  print(weighted_plot)
  grid.text("Time (days)", y = unit(0.02, "npc"), gp = gpar(fontsize = 14))

  figname_weighted = paste0('residuals-weighted', i, '.png')
  png(
    here::here('results', 'figures', figname_weighted),
    width = 8,
    height = 5,
    units = 'in',
    res = 300
  )
  print(weighted_plot)
  dev.off()
}
