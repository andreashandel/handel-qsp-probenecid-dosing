# make tables for manuscript and supplement
library(here)
library(gt)

#####################################
# best fit parameter table
#####################################
bestfit_list = readRDS(here::here('results', 'output', 'bestfit.Rds'))

nsamp = 1
# uncomment the line below to generate best fit tables for all samples of the fixed parameters
#nsamp = length(bestfit_list)
for (i in 1:nsamp) {
  bestfit = bestfit_list[[i]]
  parnames = bestfit$fitparnames
  parvalues = bestfit$solution
  parlabels = bestfit$parlabels
  tab = data.frame(parnames, parvalues, parlabels)
  # create a table using the gt package and others as needed, based on the content of the tab data frame
  gttab <- gt(tab)
  # add headers to table
  gttab <- gttab %>%
    cols_label(
      parnames = "Parameter",
      parvalues = "Value",
      parlabels = "Label"
    )

  print(gttab)

  tabname = paste0("parametertab", i, ".rds")
  # save to rds file
  saveRDS(gttab, file = here::here("results", "tables", tabname))
}
