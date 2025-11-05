##########################################
# NLL-consistent diagnostic plots (simplified)
##########################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)
library(here)
library(grid)

# ---------------- helpers ----------------

canonicalize_qty <- function(x) {
  x0 <- tolower(gsub("[^a-z0-9]", "", as.character(x)))
  dplyr::case_when(
    x0 %in%
      c("logvirusload", "viruslog10", "log10virus", "vlog10") ~ "LogVirusLoad",
    x0 %in% c("il6") ~ "IL6",
    x0 %in%
      c(
        "weight",
        "weightloss",
        "weightlossperc",
        "weightperc",
        "weightlosspercent"
      ) ~ "WeightLossPerc",
    TRUE ~ x # fallback: leave as-is
  )
}

collect_all_params <- function(bf) {
  pieces <- list(
    bf$bestpar,
    bf$bestpars,
    bf$params,
    bf$par,
    bf$fixedpars,
    if (!is.null(bf$optres)) bf$optres$solution
  )
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  unlist(pieces, use.names = TRUE)
}

gather_sigmas_tbl <- function(bestfit) {
  vec <- collect_all_params(bestfit)
  vec <- vec[!is.na(names(vec))]
  sig <- vec[grepl("^sigma_(add|prop)_", names(vec))]
  if (length(sig) == 0) {
    # no sigmas given -> empty table; downstream we treat missing as 0
    return(tibble(Quantity = character(), add = numeric(), prop = numeric()))
  }

  out <- tibble(name = names(sig), value = as.numeric(sig)) |>
    tidyr::extract(
      name,
      into = c("type", "rawqty"),
      regex = "^sigma_(add|prop)_(.+)$"
    ) |>
    mutate(Quantity = canonicalize_qty(rawqty)) |>
    select(Quantity, type, value) |>
    pivot_wider(names_from = type, values_from = value)

  # ensure columns exist even if one type was absent
  if (!"add" %in% names(out)) {
    out$add <- 0
  }
  if (!"prop" %in% names(out)) {
    out$prop <- 0
  }
  out
}

build_combined_resid_plot <- function(
  df,
  resid_col,
  ylab_txt,
  scen_labs,
  col_vals,
  shape_vals,
  var_labs
) {
  pad_df <- df %>%
    group_by(Quantity) %>%
    summarise(
      max_abs = max(abs(.data[[resid_col]]), na.rm = TRUE),
      max_abs = ifelse(is.finite(max_abs), max_abs, 0),
      ypad = max_abs * 1.1,
      x_min = min(Day, na.rm = TRUE),
      x_max = max(Day, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    transmute(Quantity, ymin = -ypad, ymax = ypad, x_min, x_max)

  p <- ggplot(
    df,
    aes(x = Day, y = .data[[resid_col]], colour = Scenario, shape = Scenario)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.9, size = 2) +
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
    ylab(ylab_txt) +
    theme_bw(base_size = 14) +
    theme(
      axis.title.x = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 12),
      legend.position = "top"
    )

  (p / plot_spacer()) + plot_layout(heights = c(1, 0.03))
}

# ---------------- main ----------------

bestfit_list <- readRDS(here::here("results", "output", "bestfit.Rds"))

nsamp <- 1
# nsamp <- length(bestfit_list)

for (i in 1:nsamp) {
  bf <- bestfit_list[[i]]
  dat <- bf$fitdata
  sim <- bf$simresult

  # ensure Day/Quantity names match model/NLL scale
  if (!("Day" %in% names(dat)) && "xvals" %in% names(dat)) {
    dat <- mutate(dat, Day = xvals)
  }
  dat <- mutate(dat, Quantity = canonicalize_qty(Quantity))

  times_to_keep <- sort(unique(dat$Day))

  # simulation on observation scale used in NLL
  sim_long <- sim %>%
    filter(time %in% times_to_keep) %>%
    select(time, Scenario, V, F, S) %>%
    mutate(V = log10(pmax(1, V))) %>%
    rename(
      Day = time,
      LogVirusLoad = V,
      IL6 = F,
      WeightLossPerc = S
    ) %>%
    pivot_longer(
      cols = c(LogVirusLoad, IL6, WeightLossPerc),
      names_to = "Quantity",
      values_to = "Predicted"
    )

  # join observed/predicted
  resid_df <- dat %>%
    inner_join(sim_long, by = c("Scenario", "Day", "Quantity")) %>%
    mutate(
      Residual = Value - Predicted,
      Quantity = factor(
        Quantity,
        levels = c("LogVirusLoad", "IL6", "WeightLossPerc")
      ),
      Scenario = factor(
        Scenario,
        levels = c("NoTreatment", "PanCytoVir10mg", "PanCytoVir100mg")
      )
    )

  # sigma table (missing types default to 0)
  sig_tbl <- gather_sigmas_tbl(bf)

  # attach sigmas, default NAs -> 0
  resid_df <- resid_df %>%
    left_join(sig_tbl, by = "Quantity") %>%
    mutate(
      add = replace_na(add, 0),
      prop = replace_na(prop, 0),
      Vi = pmax(add^2 + (prop * pmax(0, Predicted))^2, 1e-12)
    )

  # per-quantity weights (equalize across variables)
  n_by_qty <- dat %>%
    group_by(Quantity) %>%
    summarise(n = dplyr::n(), .groups = "drop")
  w_by_qty <- setNames(1 / n_by_qty$n, n_by_qty$Quantity)
  w_by_qty <- w_by_qty / mean(w_by_qty, na.rm = TRUE)

  resid_df <- resid_df %>%
    mutate(
      StdResid = Residual / sqrt(Vi),
      StdResidWeighted = sqrt(w_by_qty[as.character(Quantity)]) * StdResid
    )

  # aesthetics
  scen_labs <- c(
    NoTreatment = "No Treatment",
    PanCytoVir10mg = "10 mg/kg",
    PanCytoVir100mg = "100 mg/kg"
  )
  var_labs <- c(
    LogVirusLoad = "Log Virus Load",
    IL6 = "IL-6",
    WeightLossPerc = "Weight Loss (%)"
  )
  scen_levels <- levels(resid_df$Scenario)
  col_vals <- setNames(
    c("#0072B2", "#009E73", "#D55E00")[seq_along(scen_levels)],
    scen_levels
  ) # Okabe–Ito
  shape_vals <- setNames(c(16, 17, 15)[seq_along(scen_levels)], scen_levels)

  # ---------- (1) standardized residuals ----------
  p_std <- build_combined_resid_plot(
    resid_df,
    "StdResid",
    "Standardized Residuals",
    scen_labs,
    col_vals,
    shape_vals,
    var_labs
  )
  print(p_std)
  grid.text("Time (days)", y = unit(0.02, "npc"), gp = gpar(fontsize = 14))
  png(
    here::here(
      "results",
      "figures",
      paste0("residuals-std-combined", i, ".png")
    ),
    width = 8,
    height = 5,
    units = "in",
    res = 300
  )
  print(p_std)
  dev.off()

  # ---------- (2) weighted standardized residuals ----------
  p_stdw <- build_combined_resid_plot(
    resid_df,
    "StdResidWeighted",
    "Weighted Standardized Residuals",
    scen_labs,
    col_vals,
    shape_vals,
    var_labs
  )
  print(p_stdw)
  grid.text("Time (days)", y = unit(0.02, "npc"), gp = gpar(fontsize = 14))
  png(
    here::here(
      "results",
      "figures",
      paste0("residuals-stdweighted-combined", i, ".png")
    ),
    width = 8,
    height = 5,
    units = "in",
    res = 300
  )
  print(p_stdw)
  dev.off()

  # ---------- (2b) single-panel weighted standardized ----------
  var_colors <- c(
    LogVirusLoad = "#0072B2",
    IL6 = "#D55E00",
    WeightLossPerc = "#009E73"
  )
  df_single <- resid_df %>% select(Day, Scenario, Quantity, StdResidWeighted)
  y_max <- max(abs(df_single$StdResidWeighted), na.rm = TRUE)
  if (!is.finite(y_max) || y_max == 0) {
    y_max <- 1
  }
  y_lim <- c(-1.1 * y_max, 1.1 * y_max)

  p_single <- ggplot(
    df_single,
    aes(Day, StdResidWeighted, colour = Quantity, shape = Scenario)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.9, size = 2) +
    scale_color_manual(
      values = var_colors,
      name = "Variable",
      labels = var_labs
    ) +
    scale_shape_manual(
      values = shape_vals,
      name = "Scenario",
      labels = scen_labs
    ) +
    scale_y_continuous(limits = y_lim, expand = expansion(mult = 0)) +
    xlab("Time (days)") +
    ylab("Weighted Standardized Residuals") +
    theme_bw(base_size = 14) +
    theme(legend.position = "top")
  print(p_single)
  png(
    here::here(
      "results",
      "figures",
      paste0("residuals-stdweighted-singlepanel", i, ".png")
    ),
    width = 7.5,
    height = 5,
    units = "in",
    res = 300
  )
  print(p_single)
  dev.off()

  # ---------- (3) predicted vs observed ----------
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
    aes(Observed, Predicted, colour = Scenario, shape = Scenario)
  ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      colour = "grey50"
    ) +
    geom_point(alpha = 0.85, size = 2) +
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
  print(p_gof)
  png(
    here::here("results", "figures", paste0("pred-vs-obs-combined", i, ".png")),
    width = 8,
    height = 5,
    units = "in",
    res = 300
  )
  print(p_gof)
  dev.off()
}
