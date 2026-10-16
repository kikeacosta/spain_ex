rm (list = ls()); gc()
source("code/exploratory/00_functions_exploratory.R")
library(tidyverse)
if (!("MortalitySmooth" %in% rownames(installed.packages()))) remotes::install_github("timriffe/MortalitySmooth")
library(MortalitySmooth)

hmd <- read_rds("data_inter/deaths_exposures_hmd.rds")
cds <- unique(hmd$code) %>% sort()

for(c in cds){
  c <- "ESP"
  for(s in c("male", "female", "total")){
    plot_change(c, s, 0, 100, 1920, 2019)
    ggsave(paste0("figures/surfaces/lexis_surface_", c, "_", s, ".png"),
           w = 10, h = 10)  
  }
}

