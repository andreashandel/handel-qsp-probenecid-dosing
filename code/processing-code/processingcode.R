###############################
# processing script
#
#this script loads the raw data, processes and cleans it
#and saves it as CSV file in the processed-data folder
#
# Note the ## ---- name ---- notation
# This is done so one can pull in the chunks of code into a Quarto document
# see here: https://bookdown.org/yihui/rmarkdown-cookbook/read-chunk.html
# Not actually used here
###############################

## ---- packages --------
#load needed packages. make sure they are installed.
library(readxl) #for loading Excel files
library(dplyr) #for data processing/cleaning
library(tidyr) #for data processing/cleaning
library(here) #to set paths
library(readr) #for saving CSV

## ---- loaddata --------
#path to data
#note the use of the here() package and not absolute paths
data_location <- here::here("data", "raw-data")
virfile = paste0(data_location, "/H5N1 Plaque assay data.xlsx")
immfile = paste0(data_location, "/H5N1 ELISA serum data.xlsx")
weightfile = paste0(data_location, "/H5N1-mouse-weights.xlsx")
# load data.
virus_dat = readxl::read_excel(virfile, sheet = "DataforAnalysis")
immune_dat = readxl::read_excel(immfile, sheet = "DataforAnalysis")
weight_dat = readxl::read_excel(weightfile, sheet = "DataforAnalysis")

# only keep the scenarios we care about: uninfected/low dose/high dose
scenarios = c("NoTreatment", "PanCytoVir10mg", "PanCytoVir100mg")
dfvir <- virus_dat %>% filter(Scenario %in% scenarios)
dfimm <- immune_dat %>% filter(Scenario %in% scenarios)
dfweightchange <- weight_dat %>% filter(Scenario %in% scenarios)


## ---- combinedata --------
# virus is on natural units
# immune is on natural units
# weight is on percent change
df_all <- bind_rows(dfvir, dfimm, dfweightchange)


## ---- savedata --------
# all done, data is clean now.
# Let's assign at the end to some final variable
# makes it easier to add steps above
processeddata <- df_all
# location to save file
#save_data_location <- here::here("data","processed-data","processeddata.rds")
save_data_location <- here::here("data", "processed-data", "processeddata.csv")
readr::write_csv(processeddata, file = save_data_location)
