rm (list = ls()); gc()
source("code/exploratory/00_functions_exploratory.R")
library(tidyverse)
if (!("MortalitySmooth" %in% rownames(installed.packages()))) remotes::install_github("timriffe/MortalitySmooth")
library(MortalitySmooth)

hmd <- read_rds("data_inter/deaths_exposures_hmd.rds")
cds <- unique(hmd$code) %>% sort()

# # for(c in cds){
#   c <- "ESP"
#   for(s in c("male", "female", "total")){
#     plot_change(c, s, 0, 100, 1920, 2019)
#     # ggsave(paste0("figures/surfaces/lexis_surface_", c, "_", s, ".png"),
#     #        w = 10, h = 10)  
#   }
# # }

c <- "ESP"
s <- "female"
nplot <- "2a"
amin <- 0
amax <- 60
pmin <- 1940
pmax <- 2019

data_tst <- mort_change_data(c, s, amin, amax, pmin, pmax)

{
  data_tst
  
  qt <- 25
  
  # one for decrease of mortality (from green to blue),
  # another for mortality increases (from yellow to red)
  col_scale <- c(colorRampPalette(c("royalblue", "springgreen"), space = "Lab")(qt),
                 colorRampPalette(c("yellow", "red"), space = "Lab")(qt))
  
  # all categories
  cuts <- data_tst %>% pull(ch_cut) %>% unique() %>% sort()
  
  # assigning color to each category
  col_values <- setNames(col_scale, cuts)
  
  # Plot of mortality change over periods
  data_tst %>%
    ggplot(aes(year, age, z = ch_cut)) +
    geom_tile(aes(fill = ch_cut))+
    # adding the color palette constructed above
    scale_fill_manual(values = col_values,
                      breaks = cuts[c(1, 10, 20, 25, 30, 40, 50)],
                      name = "Mortality\nchange %")+
    #adding contour lines when the slope is 0
    geom_contour(aes(z = ch), breaks = 0, col="black",
                 alpha=.7, linewidth = .3)+
    scale_x_continuous(expand = c(0,0), breaks = seq(pmin, pmax, 10))+
    scale_y_continuous(expand = c(0,0), breaks = seq(amin, amax, 10))+
    # adding the grid of the Lexis diagram
    geom_vline(xintercept = seq(pmin, pmax, 10), linetype = "dashed",
               color = "grey30", linewidth = .20, alpha = 0.8) +
    geom_hline(yintercept = seq(amin, amax, 10), linetype = "dashed",
               color = "grey30", linewidth = .20, alpha = 0.8) +
    geom_rect(xmin = 2019.65, xmax = 2032, ymin = -1, ymax = 100.5, fill="white")+
    geom_rect(xmin = 1919, xmax = 2032, ymin = 100.05, ymax = 110.5, fill="white")+
    geom_abline(intercept = seq(-1990, -1880, 10), slope = 1,
                linetype = "dashed", color = "grey30",
                linewidth = .2, alpha = 0.8)+
    annotate("text", x = 2023, y = seq(35, 95, by = 10),
             label = as.character(rev(seq(1930, 1990, by=10))),
             size = 2, angle = 45, color="grey20") +
    annotate("text", x = seq(1980, 2020, 10), y = 103,
             label = as.character(seq(1880, 1920, by=10)),
             size = 2, angle = 45, color="grey20") +
    annotate("text", x = 2024.8, y = 28.5, label = "Cohort", size = 2.6, angle = 45, color="grey20") +
    labs(x="Period", y="Age")+
    # aesthetic details
    coord_equal(xlim = c(1920, 2032), ylim = c(-.5, 110)) +
    labs(title = paste0(c, "_", s),
         x="Period", y="Age")+
    theme_minimal()+
    # deleting default discrete legend
    theme(
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 7),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8),
      plot.title = element_text(size = 12),
      # panel.background = element_rect(fill = NA),
      # panel.grid.major = element_line(colour = NA),
      # panel.grid.minor = element_line(colour = NA),
      plot.background = element_rect(fill = "white", colour = "transparent")
    )
}
ggsave(paste0("figures/manuscript/fig0", nplot, "_lexis_surface_", c, "_", s, "_", amin, "_", amax, "_", 
              pmin, "_", pmax, ".png"),
       w = 5, h = 4)

