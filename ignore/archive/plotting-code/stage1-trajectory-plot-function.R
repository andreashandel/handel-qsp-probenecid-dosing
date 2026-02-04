################################
# Packages needed are
# ggplot2, dplyr, tidyr, patchwork, grid
################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(grid)

################################
# Stage 1 trajectory plotting function
################################
# This function builds a single multi-panel plot that overlays trajectories
# from *all* Stage 1 candidates for a small set of state variables.
#
# Expected columns in `sim_df`:
#   - time
#   - U, I, V, F, A, S (uninfected, infected, virus, innate, adaptive, morbidity)
#   - Dose
#   - Candidate
#
# Optional `data` (observations) should include:
#   - Day, Value, Quantity, Dose
#   - Quantity values should include: LogVirusLoad, IL6, WeightLossPerc
#
# The plot shows one panel per variable and overlays trajectories for every
# candidate as fine lines. Dose is mapped to color; Candidate controls grouping.
# U, I, V, and A are displayed on log10 axes, while F and S remain linear.
plot_stage1_trajectories <- function(
  sim_df,
  data = NULL,
  candidate_objectives = NULL,
  output_file = NULL,
  width = 12,
  height = 10,
  dpi = 300,
  show_plot = TRUE,
  alpha_lines = 0.25,
  data_alpha = 0.8,
  data_size = 1.5,
  data_jitter = 0.2,
  max_F = 10,
  max_S = 100,
  objective_digits = 3,
  objective_size = 2.5
) {
  if (is.null(sim_df) || !nrow(sim_df)) {
    stop("sim_df is empty; nothing to plot.")
  }

  required_cols <- c("time", "U", "I", "V", "F", "A", "S", "Dose", "Candidate")
  missing_cols <- setdiff(required_cols, names(sim_df))
  if (length(missing_cols)) {
    stop("sim_df is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Convert simulated trajectories to long format for faceting.
  plot_df <- sim_df %>%
    select(time, U, I, V, F, A, S, Dose, Candidate) %>%
    pivot_longer(
      cols = c(U, I, V, F, A, S),
      names_to = "Variable",
      values_to = "Value"
    )

  # Keep values on their original scale; log scaling is applied in panels.
  plot_df <- plot_df %>%
    mutate(
      Variable = factor(
        Variable,
        levels = c("U", "I", "V", "F", "A", "S"),
        labels = c(
          "Uninfected cells (U)",
          "Infected cells (I)",
          "Virus (V)",
          "Innate response (F)",
          "Adaptive response (A)",
          "Morbidity (S)"
        )
      ),
      Dose = factor(Dose),
      ValuePlot = Value
    )

  # Apply per-variable y-axis caps for linear panels.
  plot_df <- plot_df %>%
    mutate(
      ValuePlot = ifelse(Variable == "Innate response (F)", pmin(ValuePlot, max_F), ValuePlot),
      ValuePlot = ifelse(Variable == "Morbidity (S)", pmin(ValuePlot, max_S), ValuePlot)
    )

  # Optional data overlay for Virus, Innate response, and Morbidity.
  data_df <- NULL
  if (!is.null(data) && nrow(data)) {
    required_data_cols <- c("Day", "Value", "Quantity", "Dose")
    missing_data_cols <- setdiff(required_data_cols, names(data))
    if (length(missing_data_cols)) {
      stop("data is missing required columns: ", paste(missing_data_cols, collapse = ", "))
    }
    data_df <- data %>%
      filter(Quantity %in% c("LogVirusLoad", "IL6", "WeightLossPerc")) %>%
      mutate(
        Variable = case_when(
          Quantity == "LogVirusLoad" ~ "Virus (V)",
          Quantity == "IL6" ~ "Innate response (F)",
          Quantity == "WeightLossPerc" ~ "Morbidity (S)",
          TRUE ~ NA_character_
        ),
        ValuePlot = ifelse(Quantity == "LogVirusLoad", 10^Value, Value),
        Dose = factor(Dose)
      ) %>%
      filter(!is.na(Variable)) %>%
      mutate(
        Variable = factor(
          Variable,
          levels = levels(plot_df$Variable)
        )
      )
    # Apply the same caps to data so points stay within panel bounds.
    data_df <- data_df %>%
      mutate(
        ValuePlot = ifelse(Variable == "Innate response (F)", pmin(ValuePlot, max_F), ValuePlot),
        ValuePlot = ifelse(Variable == "Morbidity (S)", pmin(ValuePlot, max_S), ValuePlot)
      )
  }

  # Build objective labels for the morbidity panel, if provided.
  label_df <- NULL
  if (!is.null(candidate_objectives)) {
    obj_vec <- candidate_objectives
    if (is.null(names(obj_vec))) {
      # If no names, assume the vector index corresponds to candidate id.
      names(obj_vec) <- as.character(seq_along(obj_vec))
    }

    label_df <- sim_df %>%
      group_by(Candidate, Dose) %>%
      summarize(
        time = max(time),
        S = S[which.max(time)],
        .groups = "drop"
      ) %>%
      mutate(
        Objective = obj_vec[as.character(Candidate)],
        ObjectiveLabel = format(Objective, digits = objective_digits, scientific = TRUE),
        ValuePlot = pmin(S, max_S),
        Dose = factor(Dose, levels = levels(plot_df$Dose))
      ) %>%
      filter(is.finite(Objective))

    if (nrow(label_df)) {
      x_nudge <- 0.02 * max(sim_df$time)
      label_df$time <- label_df$time + x_nudge
    }
  }

  # Build a single-variable panel so each variable can have its own scale.
  build_panel <- function(df_line,
                          df_point = NULL,
                          log_scale = FALSE,
                          y_max = NULL,
                          label_df = NULL) {
    p <- ggplot(
      df_line,
      aes(
        x = time,
        y = ValuePlot,
        colour = Dose,
        group = interaction(Dose, Candidate)
      )
    ) +
      geom_line(linewidth = 0.5, alpha = alpha_lines) +
      {
        if (!is.null(df_point) && nrow(df_point)) {
          geom_point(
            data = df_point,
            aes(x = Day, y = ValuePlot, colour = Dose),
            size = data_size,
            alpha = data_alpha,
            position = position_jitter(width = data_jitter, height = 0),
            inherit.aes = FALSE
          )
        } else {
          NULL
        }
      } +
      {
        if (!is.null(label_df) && nrow(label_df)) {
          geom_text(
            data = label_df,
            aes(x = time, y = ValuePlot, label = ObjectiveLabel, colour = Dose),
            size = objective_size,
            hjust = 0,
            inherit.aes = FALSE
          )
        } else {
          NULL
        }
      } +
      scale_x_continuous(expand = expansion(mult = c(0.02, 0.12))) +
      labs(
        x = "Time (days)",
        y = "Value",
        colour = "Dose"
      ) +
      theme_bw() +
      theme(
        aspect.ratio = 1,
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(size = 10)
      )

    if (log_scale) {
      p <- p + scale_y_log10(
        limits = c(1e-5, NA),
        oob = scales::oob_squish,
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      )
    }
    if (!is.null(y_max)) {
      p <- p + coord_cartesian(ylim = c(0, y_max))
    }

    p
  }

  make_var_panel <- function(var_label, log_scale = FALSE, y_max = NULL, label_df = NULL) {
    df_line <- plot_df %>% filter(Variable == var_label)
    df_point <- if (!is.null(data_df)) data_df %>% filter(Variable == var_label) else NULL
    build_panel(df_line, df_point, log_scale = log_scale, y_max = y_max, label_df = label_df) +
      ggtitle(var_label)
  }

  p_U <- make_var_panel("Uninfected cells (U)", log_scale = TRUE)
  p_I <- make_var_panel("Infected cells (I)", log_scale = TRUE)
  p_V <- make_var_panel("Virus (V)", log_scale = TRUE)
  p_F <- make_var_panel("Innate response (F)", log_scale = FALSE, y_max = max_F)
  p_A <- make_var_panel("Adaptive response (A)", log_scale = TRUE)
  p_S <- make_var_panel("Morbidity (S)", log_scale = FALSE, y_max = max_S, label_df = label_df)

  # Arrange panels: top row U, I, V; bottom row F, A, S.
  p <- patchwork::wrap_plots(list(p_U, p_I, p_V, p_F, p_A, p_S), ncol = 3)

  # Save to file if requested
  if (!is.null(output_file)) {
    ggsave(
      filename = output_file,
      plot = p,
      width = width,
      height = height,
      dpi = dpi
    )
  }

  # Display plot interactively if requested
  if (show_plot) {
    if (!interactive()) {
      grDevices::dev.new()
    }
    print(p)
    grDevices::dev.flush()
  }

  invisible(p)
}
