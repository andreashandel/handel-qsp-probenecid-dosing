################################
# Packages needed are
# ggplot2, dplyr, patchwork
################################
library(ggplot2)
library(dplyr)
library(patchwork)

################################
# Time-series plotting function
################################
plot_timeseries <- function(
  data = NULL,
  modelfit,
  tmax = 7,
  dose_levels = NULL,
  dose_levels_labels = NULL,
  x_jitter = 0,
  # Optional RNG seed so jitter is reproducible; set to a number (e.g., 123) or leave NULL.
  x_jitter_seed = 1234
) {
  ## -------------------------------------------------------------------------
  ## Aesthetics and labels
  ## -------------------------------------------------------------------------
  okabe_ito <- c(
    "#0072B2",
    "#009E73",
    "#D55E00",
    "#E69F00",
    "#56B4E9",
    "#F0E442",
    "#CC79A7",
    "#000000"
  )
  color_vals <- NULL
  linetype_base <- c("solid", "42", "12", "13", "F4", "F2")
  shape_base <- c(16, 17, 15, 18, 19, 8)

  # --- NEW: build a reusable position object for points (x-only jitter) ------
  # If x_jitter == 0, no jitter is applied. Otherwise, apply horizontal jitter
  # with the requested width and no vertical jitter. Using 'seed' makes jitter
  # reproducible across calls/panels.
  point_pos <- if (x_jitter > 0) {
    position_jitter(width = x_jitter, height = 0, seed = x_jitter_seed)
  } else {
    position_identity()
  }

  ## determine dose ordering and labels --------------------------------------
  infer_doses <- function(df) {
    if (!is.null(df) && nrow(df) > 0 && "Dose" %in% names(df)) {
      return(df$Dose)
    }
    numeric(0)
  }

  if (is.null(dose_levels)) {
    candidate_levels <- sort(unique(c(
      infer_doses(modelfit),
      infer_doses(data)
    )))
    if (!length(candidate_levels)) {
      stop(
        "Unable to infer dose levels; please provide `dose_levels` explicitly."
      )
    }
    dose_levels <- candidate_levels
  }

  dose_levels <- sort(unique(as.numeric(dose_levels)))

  if (is.null(dose_levels_labels)) {
    dose_levels_labels <- vapply(
      dose_levels,
      function(d) {
        if (isTRUE(all.equal(d, 0))) {
          "no drug"
        } else {
          paste0(format(d, trim = TRUE, scientific = FALSE), " mg/kg")
        }
      },
      character(1)
    )
  } else if (length(dose_levels_labels) != length(dose_levels)) {
    stop("`dose_levels` and `dose_levels_labels` must have the same length.")
  }

  color_vals <- rep(okabe_ito, length.out = length(dose_levels_labels))
  linetype_vals <- rep(linetype_base, length.out = length(dose_levels_labels))
  shape_vals <- rep(shape_base, length.out = length(dose_levels_labels))

  recode_dose <- function(df) {
    if (is.null(df) || !nrow(df)) {
      return(df)
    }
    df %>%
      mutate(
        Dose = factor(Dose, levels = dose_levels, labels = dose_levels_labels)
      )
  }

  ## -------------------------------------------------------------------------
  ## Ensure Dose is a factor with the desired labels
  ## -------------------------------------------------------------------------
  modelfit <- recode_dose(modelfit)

  if (!is.null(data) && nrow(data) > 0) {
    data <- recode_dose(data)

    Vvals <- data %>%
      filter(Quantity == "LogVirusLoad") %>%
      mutate(Value = 10^Value)
    Innvals <- data %>% filter(Quantity == "IL6")
    Symvals <- data %>% filter(Quantity == "WeightLossPerc")
  } else {
    Vvals <- Innvals <- Symvals <- NULL
  }

  ## -------------------------------------------------------------------------
  ## Determine dynamic y-axis limits so all trajectories remain visible
  ## -------------------------------------------------------------------------
  collect_values <- function(var_name, df_point = NULL) {
    line_vals <- modelfit[[var_name]]
    extra_vals <- if (!is.null(df_point)) df_point$Value else numeric(0)
    vals <- c(line_vals, extra_vals)
    vals[is.finite(vals)]
  }

  compute_limits <- function(
    var_name,
    logy = FALSE,
    df_point = NULL,
    clamp_zero = FALSE
  ) {
    vals <- collect_values(var_name, df_point)
    if (!length(vals)) {
      return(NULL)
    }

    if (logy) {
      pos_vals <- vals[vals > 0]
      if (!length(pos_vals)) {
        warning(
          sprintf(
            "All values for '%s' are non-positive; using a default log10 range.",
            var_name
          )
        )
        return(c(1e-8, 1))
      }
      min_val <- min(pos_vals)
      max_val <- max(pos_vals)
      if (min_val == max_val) {
        lower <- min_val / 2
        upper <- max_val * 2
      } else {
        lower <- min_val / 1.5
        upper <- max_val * 1.25
      }
      lower <- max(lower, min_val / 10)
      lower <- max(lower, 1e-8)
      upper <- max(upper, lower * 1.1)
      c(lower, upper)
    } else {
      min_val <- min(vals)
      max_val <- max(vals)
      if (min_val == max_val) {
        pad <- ifelse(min_val == 0, 1, abs(min_val) * 0.1)
      } else {
        pad <- (max_val - min_val) * 0.1
      }
      if (!is.finite(pad) || pad == 0) {
        pad <- 1
      }
      lower <- min_val - pad
      upper <- max_val + pad
      if (clamp_zero && min_val >= 0) {
        lower <- 0
      }
      if (lower == upper) {
        lower <- lower - 1
        upper <- upper + 1
      }
      c(lower, upper)
    }
  }

  uc_limits <- compute_limits("U", logy = TRUE)
  ic_limits <- compute_limits("I", logy = TRUE)
  vir_limits <- compute_limits("V", logy = TRUE, df_point = Vvals)
  inn_limits <- compute_limits("F", df_point = Innvals, clamp_zero = TRUE)
  ada_limits <- compute_limits("A", logy = TRUE)
  sym_limits <- compute_limits("S", df_point = Symvals, clamp_zero = TRUE)
  ad_limits <- compute_limits("Ad", logy = TRUE)
  ac_limits <- compute_limits("Ac", logy = TRUE)
  at_limits <- compute_limits("At", logy = TRUE)

  ## -------------------------------------------------------------------------
  ## Re-usable panel builder
  ## -------------------------------------------------------------------------
  plot_template <- function(
    df_line,
    df_point = NULL,
    xvar,
    yvar,
    ylabel,
    logy = FALSE,
    ylimits = NULL,
    tmax = 7,
    keep_legend = FALSE
  ) {
    # define the guide only when this panel should provide it
    colour_guide <- if (keep_legend) {
      guide_legend(
        override.aes = list(
          linetype = linetype_vals,
          shape = shape_vals,
          size = 1.1
        )
      )
    } else {
      "none"
    }

    p <- ggplot() +
      geom_line(
        data = df_line,
        aes(x = !!sym(xvar), y = !!sym(yvar), colour = Dose, linetype = Dose),
        linewidth = 1.2,
        lineend = "round",
        show.legend = keep_legend # <- only this panel feeds the legend
      ) +
      labs(x = NULL, y = ylabel) +

      ## ---- single combined legend comes from the colour scale --------------
      scale_colour_manual(
        name = "Dose", # legend title
        values = setNames(color_vals, dose_levels_labels),
        breaks = dose_levels_labels,
        guide = colour_guide # <- legend only on the chosen panel
      ) +
      ## hide separate linetype & shape legends
      scale_linetype_manual(
        values = setNames(linetype_vals, dose_levels_labels),
        guide = "none"
      ) +
      scale_shape_manual(
        values = setNames(shape_vals, dose_levels_labels),
        guide = "none"
      )

    ## add points if measurement data exist
    if (!is.null(df_point) && nrow(df_point) > 0) {
      p <- p +
        geom_point(
          data = df_point,
          aes(x = xvals, y = Value, colour = Dose, shape = Dose),
          size = 1,
          alpha = 0.5,
          show.legend = FALSE, # no separate shape legend
          position = point_pos # <-- NEW: x-only jitter applied here
        )
    }

    ## y-axis scaling
    if (logy) {
      if (is.null(ylimits)) {
        p <- p + scale_y_log10(expand = expansion(mult = c(0, 0.05)))
      } else {
        p <- p +
          scale_y_log10(limits = ylimits, expand = expansion(mult = c(0, 0.05)))
      }
    } else if (!is.null(ylimits)) {
      p <- p +
        scale_y_continuous(
          limits = ylimits,
          expand = expansion(mult = c(0, 0.05))
        )
    }

    ## x-axis settings
    p <- p +
      scale_x_continuous(
        limits = c(0, tmax),
        breaks = seq(0, tmax, 1),
        minor_breaks = NULL
      ) +
      theme_minimal() +
      theme(
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.key.width = unit(2, "cm"),
        axis.title.y = element_text(size = 10),
        legend.position = "none" # legends collected by patchwork
      )

    return(p)
  }

  ## -------------------------------------------------------------------------
  ## Build the panels, one for each model variable
  ## -------------------------------------------------------------------------
  p_uc <- plot_template(
    modelfit,
    NULL,
    "time",
    "U",
    "Uninfected Cells",
    logy = TRUE,
    ylimits = uc_limits,
    tmax = tmax,
    keep_legend = TRUE
  )
  p_ic <- plot_template(
    modelfit,
    NULL,
    "time",
    "I",
    "Infected Cells",
    logy = TRUE,
    ylimits = ic_limits,
    tmax = tmax
  )
  p_vir <- plot_template(
    modelfit,
    Vvals,
    "time",
    "V",
    "Virus Load",
    logy = TRUE,
    #ylimits = c(1e-2, vir_limits[2]),
    ylimits = vir_limits,
    tmax = tmax
  )
  p_inn <- plot_template(
    modelfit,
    Innvals,
    "time",
    "F",
    "Innate Response (IL-6)",
    logy = FALSE,
    ylimits = inn_limits,
    tmax = tmax
  )
  p_ada <- plot_template(
    modelfit,
    NULL,
    "time",
    "A",
    "Adaptive Response",
    logy = TRUE,
    ylimits = ada_limits,
    tmax = tmax
  )
  p_sym <- plot_template(
    modelfit,
    Symvals,
    "time",
    "S",
    "Morbidity (weight loss)",
    logy = FALSE,
    ylimits = sym_limits,
    tmax = tmax
  )
  p_ad <- plot_template(
    modelfit,
    NULL,
    "time",
    "Ad",
    "Drug depot",
    logy = TRUE,
    ylimits = ad_limits,
    tmax = tmax
  )
  p_ac <- plot_template(
    modelfit,
    NULL,
    "time",
    "Ac",
    "Drug central",
    logy = TRUE,
    ylimits = ac_limits,
    tmax = tmax
  )
  p_at <- plot_template(
    modelfit,
    NULL,
    "time",
    "At",
    "Drug target",
    logy = TRUE,
    ylimits = at_limits,
    tmax = tmax
  )

  ## -------------------------------------------------------------------------
  ## Assemble, collect legend, add global x-axis label
  ## -------------------------------------------------------------------------
  combined_plot <- (p_uc | p_ic | p_vir) /
    (p_inn | p_ada | p_sym) /
    (p_ad | p_ac | p_at) +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "center"
    )

  final_plot <- combined_plot +
    plot_annotation(
      caption = "Days",
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

  return(final_plot)
}
