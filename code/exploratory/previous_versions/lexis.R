rm(list=ls())
# source("code/00_functions.R")
library(tidyverse)

amin <- 0
amax <- 110
pmin <- 1930
pmax <- 2030
bt_coh <- 1940

pol1 <- c(0, 10, 10, 0)
pol2 <- c(0, 10, 10)
pol3 <- c(0, 10, 0)
ags <- seq(0, 90, 10)

polin1 <- polin2 <- polin3 <- NULL
for(i in seq(0, 90, 10)) polin1 <- c(polin1, pol1+i)
for(i in seq(0, 90, 10)) polin2 <- c(polin2, pol2+i)
for(i in seq(0, 90, 10)) polin3 <- c(polin3, pol3+i)

parallelogram_data <- data.frame(
  x = rep(c(1940, 1950, 2020, 2010), 10),
  y = polin1,
  id = rep(1:10, each = 4)
)

parallelogram_data2 <- data.frame(
  x = rep(c(1940, 1950, 1940), 10),
  y = polin2,
  id = rep(1:10, each = 3)
)

parallelogram_data3 <- data.frame(
  x = rep(c(2010, 2020, 2020), 10),
  y = polin3,
  id = rep(11:20, each = 3)
)

parallelogram_data4 <- 
  bind_rows(parallelogram_data2,
            parallelogram_data3)

ggplot()+
  # adding lines for periods
  geom_vline(xintercept = seq(pmin, pmax, 10), 
             linewidth = 0.2, linetype = "dashed", alpha = 0.8, color = "grey30")+
  # adding lines for ages
  geom_hline(yintercept = seq(amin, amax, 10), 
             linewidth = 0.2, linetype = "dashed", alpha = 0.8, color = "grey30")+
  geom_abline(intercept = seq(-pmax, -(pmin-amax), 10), slope = 1, 
              linetype = "dashed", color = "grey30", linewidth = .2, alpha = 0.8)+
  # geom_vline(xintercept = c(1940, 2020))+
  geom_polygon(
    data = parallelogram_data,  # Pass the separate data frame
    aes(x = x, y = y, group = id),         # Map the coordinates
    fill  = "blue",             # Fill color
    alpha = 0.2,               # Transparency
    color = "black",           # Optional: Border color
    linewidth = .3              # Optional: Border thickness
  )+
  geom_polygon(
    data = parallelogram_data4,  # Pass the separate data frame
    aes(x = x, y = y, group = id),         # Map the coordinates
    fill  = "red",             # Fill color
    alpha = 0.2,               # Transparency
    color = "black",           # Optional: Border color
    linewidth = .3              # Optional: Border thickness
  )+
  # geom_text(aes(x = 2022, y = 85, label = "1940"), angle = 45, hjust = 0, size = 3)+
  geom_rect(aes(xmin = 2030.5, xmax = 2040, ymin = -.5, ymax = 120), fill = "white")+
  # geom_text(aes(x = 2022, y = seq(25, 95, 10), label = rev(seq(1930, 2000, 10))), angle = 45, hjust = 0, size = 2.5)+
  geom_text(aes(x = 2032, y = seq(22, 102, 10), label = rev(seq(1930, 2010, 10))), 
            angle = 45, hjust = 0, size = 3)+
  geom_text(aes(x = 2031, y = 15, label = "Cohort"), 
            angle = 45, hjust = 0, size = 3.5)+
  scale_x_continuous(breaks = seq(1910, 2020, 10))+
  scale_y_continuous(breaks = seq(0, 100, 10))+
  coord_fixed(expand = 0, xlim = c(1930, 2040), ylim = c(0, 110))+
  labs(y = "Age", x = "Period")+
  theme_minimal()+
  theme(
    # panel.border = element_blank(),
    #     plot.background = element_rect(fill = "white")
    )

ggsave("figures/lexis_tst2.png", dpi = 600,
       w = 5, h = 5)

