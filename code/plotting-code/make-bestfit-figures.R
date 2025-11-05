##########################################
# make best fit time-series plot
# makes one figure for each sample of the fixed parameters
##########################################
# packages needed by this script and functions it calls
library(ggplot2)
library(dplyr)
library(patchwork)
library(here)

source(here::here('code/plotting-code/timeseries-plot-function.R'))
bestfit_list = readRDS(here::here('results', 'output', 'bestfit.Rds'))
#bestfit_list = readRDS(here::here('results', 'output', 'bestfit-old.Rds'))

dose_levels_labels <- c("no drug", "10 mg/kg", "100 mg/kg")
dose_levels <- c(0, 10, 100)

nsamp = 1
# uncomment the line below to generate best fit figures for all samples of the fixed parameters
#nsamp = length(bestfit_list)
for (i in 1:nsamp) {
  bestfit <- bestfit_list[[i]]
  bestfit_plots <- plot_timeseries(
    data = bestfit$fitdata,
    modelfit = bestfit$simresult,
    tmax = 7,
    dose_levels,
    dose_levels_labels,
    x_jitter = 0.3
  )
  plot(bestfit_plots)
  # save plot as png file into the results folder
  figname = paste0('bestfit', i, '.png')
  png(
    here::here('results', 'figures', figname),
    width = 8,
    height = 5,
    units = 'in',
    res = 300
  )
  plot(bestfit_plots)
  dev.off()
}
