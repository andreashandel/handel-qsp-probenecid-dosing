##############################
# function to plot dose-response figures from a list of data frames
##############################
# packages needed by this function are
# ggplot2, dplyr, patchwork, tidyr
# those are loaded in the calling scripts

plot_outcomes <- function(df_list, scenarios) {
  ## --- Colour-blind-safe palette (blue, bluish-green, vermillion) ----------
  base_pal <- c("#0072B2", "#009E73", "#D55E00")

  # keep only the 'reduction_df' from each element
  df_list <- lapply(df_list, `[[`, "reduction_df")

  ## Combine list -> single long data frame
  df_all <- bind_rows(df_list, .id = "rep") %>% # 'rep' not used, but keeps provenance
    filter(Scenario %in% scenarios) %>%
    mutate(Scenario = factor(Scenario, levels = scenarios))

  ## Long format for generic summarization across outcome variables
  long <- df_all %>%
    pivot_longer(
      cols = c(perc_AUCV, perc_AUCF, perc_AUCS),
      names_to = "metric",
      values_to = "value"
    )

  ## get baseline from the first list entry and reshape like `long`
  baseline_long <- df_list[[1]] %>%
    dplyr::ungroup() %>%
    filter(Scenario %in% scenarios) %>%
    mutate(Scenario = factor(Scenario, levels = scenarios)) %>%
    pivot_longer(
      cols = c(perc_AUCV, perc_AUCF, perc_AUCS),
      names_to = "metric",
      values_to = "baseline"
    ) %>%
    select(Scenario, Dose, metric, baseline)

  ## Summary stats across replicates in df_list
  ## 95% bounds as 2.5% and 97.5% quantiles
  ## also attach results for baseline values of fixed parameters
  summ <- long %>%
    group_by(Scenario, Dose, metric) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      lower = quantile(value, 0.025, na.rm = TRUE, type = 8),
      upper = quantile(value, 0.975, na.rm = TRUE, type = 8),
      .groups = "drop"
    ) %>%
    left_join(baseline_long, by = c("Scenario", "Dose", "metric"))

  ## Named vectors → each Scenario gets its own hue/linetype
  col_vals <- setNames(base_pal[seq_along(scenarios)], scenarios)
  linetype_vals <- setNames(
    rep(c("solid", "dashed", "dotted"), length.out = length(scenarios)),
    scenarios
  )

  axis_title_size <- 14
  axis_text_size <- 12

  # Base theme
  base_theme <- theme_minimal() +
    theme(
      plot.title = element_blank(),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.key.width = unit(2, "cm"),
      axis.title = element_text(size = axis_title_size),
      axis.text = element_text(size = axis_text_size),
      plot.margin = margin(t = 5, r = 15, b = 5, l = 5)
    )

  # Helper to build each panel from summarised data
  make_plot <- function(metric_name, ylab, keep_legend = FALSE) {
    dfm <- filter(summ, metric == metric_name)

    p <- ggplot(
      dfm,
      aes(x = Dose, y = baseline, colour = Scenario, linetype = Scenario)
    ) +
      geom_ribbon(
        aes(ymin = lower, ymax = upper, fill = Scenario),
        alpha = 0.2,
        colour = NA,
        show.legend = FALSE
      ) +
      geom_line(linewidth = 1.2) +
      # use below if I also want to plot the mean in addition to the baseline
      # geom_line(
      #   aes(y = mean),
      #   linewidth = 0.9,
      #   linetype = "dotdash",
      #   show.legend = FALSE,
      #   na.rm = TRUE
      # ) +
      scale_color_manual(values = col_vals) +
      scale_fill_manual(values = col_vals) +
      scale_linetype_manual(values = linetype_vals) +
      scale_x_log10(breaks = c(1e-3, 1e-1, 1e1, 1e3)) +
      labs(
        x = NULL,
        y = ylab,
        colour = "Scenario:",
        linetype = "Scenario:"
      ) +
      base_theme +
      geom_vline(
        xintercept = c(10, 100),
        linetype = "dashed",
        colour = "black"
      )

    if (!keep_legend) {
      p <- p + theme(legend.position = "none")
    }
    p
  }

  ## Only the first panel keeps its legend; patchwork will collect it
  p1 <- make_plot(
    "perc_AUCV",
    "Log Viral Load Reduction (%)",
    keep_legend = TRUE
  )
  p2 <- make_plot("perc_AUCF", "Innate Response Reduction (%)")
  p3 <- make_plot("perc_AUCS", "Morbidity Reduction (%)")

  ## --- Combine, collect guides, place legend on top & add global x-label ---
  combined <- (p1 | p2 | p3) +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center"
    )

  combined +
    plot_annotation(
      caption = "Dose",
      theme = theme(
        plot.caption = element_text(
          hjust = 0.5,
          vjust = -0.5,
          face = "bold",
          size = 14
        ),
        plot.margin = margin(b = 20)
      )
    )
}
