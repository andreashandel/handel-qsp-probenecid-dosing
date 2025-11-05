################################
# Visual predictive check plotting utilities
################################
# Required packages: ggplot2, dplyr, tidyr

library(ggplot2)
library(dplyr)
library(tidyr)

#' Build a visual predictive check (VPC) plot for a given dosing schedule.
#'
#' @param sim_timeseries List of data frames, each containing simulated
#'   time-series output (one element per stochastic replicate). Every data frame
#'   must contain the columns `Schedule`, `Dose`, `time`, `V`, `F`, and `S`.
#' @param schedule_id Character identifier of the schedule to visualise.
#' @param observed_df Data frame with observed measurements. Must include the
#'   columns `Dose`, `Quantity`, `Day`, and `Value`.
#' @param dose_levels Optional numeric vector defining the ordering of dose
#'   levels. If `NULL`, the levels are inferred from the simulation results.
#' @param dose_labels Optional character vector (same length as `dose_levels`)
#'   that provides pretty labels for the legend.
#' @param quantiles Numeric vector of length three containing the lower,
#'   median, and upper quantiles to display.
#' @return A `ggplot2` object representing the VPC for the requested schedule.
plot_vpc_timeseries <- function(
  sim_timeseries,
  schedule_id,
  observed_df,
  dose_levels = NULL,
  dose_labels = NULL,
  quantiles = c(0.05, 0.5, 0.95)
) {
  if (length(quantiles) != 3L) {
    stop("`quantiles` must contain exactly three probabilities (lower, median, upper).")
  }

  lower_q <- quantiles[1]
  median_q <- quantiles[2]
  upper_q <- quantiles[3]

  if (!lower_q < median_q || !median_q < upper_q) {
    stop("`quantiles` must be strictly increasing (e.g., c(0.05, 0.5, 0.95)).")
  }

  sim_schedule <- lapply(sim_timeseries, function(df) {
    df %>%
      filter(Schedule == schedule_id)
  })

  sim_schedule <- sim_schedule[vapply(sim_schedule, nrow, integer(1)) > 0]

  if (!length(sim_schedule)) {
    stop(sprintf("No simulated trajectories available for schedule '%s'.", schedule_id))
  }

  sim_long <- bind_rows(sim_schedule, .id = "rep") %>%
    select(rep, Schedule, Dose, time, V, F, S) %>%
    mutate(
      Dose = suppressWarnings(as.numeric(as.character(Dose)))
    ) %>%
    pivot_longer(
      cols = c(V, F, S),
      names_to = "Variable",
      values_to = "SimValue"
    ) %>%
    mutate(
      Quantity = recode(Variable, V = "LogVirusLoad", F = "IL6", S = "Weight"),
      Value = if_else(Quantity == "LogVirusLoad", log10(pmax(SimValue, 1)), SimValue)
    )

  if (is.null(dose_levels)) {
    dose_levels <- sort(unique(sim_long$Dose))
  }

  if (is.null(dose_labels)) {
    format_label <- function(d) {
      if (isTRUE(all.equal(d, 0))) {
        "no drug"
      } else {
        paste0(format(d, trim = TRUE, scientific = FALSE), " mg/kg")
      }
    }
    dose_labels <- vapply(dose_levels, format_label, character(1))
  }

  dose_palette <- c("#0072B2", "#009E73", "#D55E00", "#E69F00", "#56B4E9")
  dose_palette <- rep(dose_palette, length.out = length(dose_labels))
  names(dose_palette) <- dose_labels

  sim_summary <- sim_long %>%
    filter(!is.na(Dose)) %>%
    mutate(
      Dose = factor(Dose, levels = dose_levels, labels = dose_labels)
    ) %>%
    group_by(Dose, Quantity, time) %>%
    summarise(
      lower = quantile(Value, lower_q, na.rm = TRUE, type = 8),
      median = quantile(Value, median_q, na.rm = TRUE, type = 8),
      upper = quantile(Value, upper_q, na.rm = TRUE, type = 8),
      .groups = "drop"
    )

  obs_schedule <- observed_df %>%
    filter(Dose %in% dose_levels) %>%
    mutate(
      Dose = factor(Dose, levels = dose_levels, labels = dose_labels)
    )

  var_labels <- c(
    LogVirusLoad = "Log Virus Load",
    IL6 = "IL-6",
    Weight = "Weight Loss"
  )

  ggplot(sim_summary, aes(x = time, colour = Dose, fill = Dose)) +
    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      colour = NA,
      alpha = 0.2
    ) +
    geom_line(aes(y = median), linewidth = 1) +
    geom_point(
      data = obs_schedule,
      aes(x = Day, y = Value),
      size = 2,
      alpha = 0.8,
      inherit.aes = FALSE
    ) +
    facet_wrap(
      ~Quantity,
      scales = "free_y",
      labeller = labeller(Quantity = var_labels)
    ) +
    scale_colour_manual(values = dose_palette, name = "Dose") +
    scale_fill_manual(values = dose_palette, guide = "none") +
    labs(x = "Time (days)", y = NULL) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "top",
      strip.background = element_blank(),
      strip.text = element_text(size = 12)
    )
}
