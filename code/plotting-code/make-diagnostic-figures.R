##########################################
# make residual plot for best fit  
##########################################
############################################
##  Packages
############################################
library(tidyverse)   # dplyr, tidyr, ggplot2, readr, …
library(patchwork)   # for combining ggplot objects
library(here)       # for file paths
library(grid)        # for a custom x-axis label


############################################
##  1.  Read in the fit results
############################################

bestfit_list = readRDS(here::here('results','output','bestfit.Rds'))


nsamp = 1
# uncomment the line below to generate best fit tables for all samples of the fixed parameters 
nsamp = length(bestfit_list)
for (i in 1:nsamp) 
{
  bestfit = bestfit_list[[i]]
  dat  <- bestfit$fitdata
  sim  <- bestfit$simresult


times_to_keep <- unique(dat$xvals)

sim_long <- sim %>%
  filter(time %in% times_to_keep) %>%      # keep only matching times
  select(time, Scenario, V, F, S) %>%  # keep relevant columns
  mutate(
    V =log10(V + 1e-12)
  ) %>%
  rename(Day = time,
        LogVirusLoad = V,
         IL6          = F,
         Weight       = S) %>%
  pivot_longer(
    cols      = c(LogVirusLoad, IL6, Weight),
    names_to  = "Quantity",
    values_to = "Predicted"
  )


###################################
# 3.  One-to-one join & residuals
###################################
resid_df <- dat %>%                       # keep every data row (incl. replicates)
  inner_join(sim_long,                   # exact match on variable, time, scenario
             by = c("Scenario", "Day", "Quantity")) %>%
  mutate(
    Residual = Value - Predicted,
    # facet order
    Quantity = factor(Quantity,
                      levels = c("LogVirusLoad", "IL6", "Weight")),
    Scenario = factor(Scenario,
                      levels = c("NoTreatment",
                                 "PanCytoVir10mg",
                                 "PanCytoVir100mg"))
  )



###################################
# 4.  Nice facet labels
###################################
scen_labs <- c(
  NoTreatment       = "No Treatment",
  PanCytoVir10mg    = "10 mg/kg",
  PanCytoVir100mg   = "100 mg/kg"
)
var_labs  <- c(
  LogVirusLoad = "Log Virus Load",
  IL6          = "IL-6",
  Weight       = "Weight"
)

###################################
# 5.  Build the 3×3 residual plot
###################################
p <- ggplot(resid_df,
            aes(x = Day, y = Residual)) +
  geom_hline(yintercept = 0,
             linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.65, size = 2) +
  facet_grid(
    rows  = vars(Quantity),
    cols  = vars(Scenario),
    labeller = labeller(Quantity = var_labs,
                        Scenario = scen_labs),
    scales = "free_y"          # tighter scale for IL-6 row
  ) +
  ylab("Residual") +
  theme_bw(base_size = 14) +
  theme(
    axis.title.x = element_blank(),   # master label added later
    strip.background = element_blank(),
    strip.text = element_text(size = 12)
  )

###################################
# 6.  Add a single centred x-axis
###################################
#   Patchwork trick: stack a blank
#   spacer row, then draw the label
###################################
blank <- plot_spacer()

full_plot <-
  p / blank +                            # stack blank row below
  plot_layout(heights = c(1, 0.03)) &    # tiny height for spacer
  theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5))

# draw the plot and then the label
print(full_plot)
grid.text("Time (days)",
          y = unit(0.02, "npc"),
          gp = gpar(fontsize = 14))

# save to png file
  figname = paste0('residuals',i,'.png')
  png(here::here('results','figures',figname), width = 8, height = 5, units = 'in', res = 300)
  print(full_plot)
  dev.off()
}

