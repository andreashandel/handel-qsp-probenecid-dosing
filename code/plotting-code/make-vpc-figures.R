##########################################
# Create visual predictive check (VPC) plots
##########################################

library(dplyr)
library(here)

source(here::here("code", "plotting-code", "vpc-plot-function.R"))

# Load simulation outputs (list of stochastic replicates)
dose_response_list <- readRDS(here::here("results", "output", "dose-response-results.Rds"))

# Extract the time-series component for every replicate
sim_timeseries <- lapply(dose_response_list, `[[`, "timeseries_df")
if (!length(sim_timeseries)) {
  stop("No time-series trajectories found in dose-response results.")
}

# Determine available schedules from the first replicate that has data
first_with_data <- which(vapply(sim_timeseries, function(df) nrow(df) > 0, logical(1)))[1]
if (is.na(first_with_data)) {
  stop("No schedule information found in the simulation outputs.")
}

schedule_ids <- sort(unique(sim_timeseries[[first_with_data]]$Schedule))

# Observed data for overlay: take the first best-fit sample
bestfit_list <- readRDS(here::here("results", "output", "bestfit.Rds"))
if (!length(bestfit_list)) {
  stop("Best-fit object not found or empty; cannot build VPC plots.")
}

observed_data <- bestfit_list[[1]]$fitdata %>%
  mutate(
    Day = xvals,
    Dose = dplyr::case_when(
      Scenario == "NoTreatment" ~ 0,
      Scenario == "PanCytoVir10mg" ~ 10,
      Scenario == "PanCytoVir100mg" ~ 100,
      TRUE ~ NA_real_ # NOTE: update this mapping if additional scenarios are introduced.
    )
  ) %>%
  filter(!is.na(Dose)) %>%
  filter(Quantity %in% c("LogVirusLoad", "IL6", "Weight")) %>%
  select(Dose, Quantity, Day, Value)

if (!nrow(observed_data)) {
  warning("No observed data matched the known dose levels; VPC plots will omit observations.")
}

# Pretty labels for the legend (consistent with other figures)
format_dose_label <- function(dose) {
  if (isTRUE(all.equal(dose, 0))) {
    return("no drug")
  }
  paste0(format(dose, trim = TRUE, scientific = FALSE), " mg/kg")
}

all_dose_levels <- observed_data$Dose %>%
  unique() %>%
  sort()

all_dose_labels <- vapply(all_dose_levels, format_dose_label, character(1))
if (!length(all_dose_levels)) {
  all_dose_levels <- NULL
  all_dose_labels <- NULL
}

# Generate and save VPC plots for each schedule
for (schedule_id in schedule_ids) {
  vpc_plot <- plot_vpc_timeseries(
    sim_timeseries = sim_timeseries,
    schedule_id = schedule_id,
    observed_df = observed_data,
    dose_levels = all_dose_levels,
    dose_labels = all_dose_labels
  )

  print(vpc_plot)

  outfile <- here::here(
    "results",
    "figures",
    sprintf("vpc-%s.png", schedule_id)
  )

  ggsave(outfile, vpc_plot, width = 10, height = 6, dpi = 300)
}
