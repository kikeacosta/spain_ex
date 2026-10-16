# ==============================================================================
# 04. Figure 2. Lexis diagram of the panel design used to decompose life
#     expectancy convergence into cohort-specific contributions (1940-2019,
#     10-year age groups). Blue parallelograms: cohorts observed at every age
#     of the group (kept); red triangles: cohorts partially observed (dropped).
# Output: figures/manuscript/fig02_lexis_diagram.png
# ==============================================================================
rm(list = ls()); gc()
library(tidyverse)

amin <- 0;    amax <- 110
pmin <- 1930; pmax <- 2030
ymin <- 1940; ymax <- 2020 # observation window

# vertices of the parallelograms (complete cohorts) and triangles (partial)
pol1 <- c(0, 10, 10, 0)
pol2 <- c(0, 10, 10)
pol3 <- c(0, 10, 0)

polin1 <- polin2 <- polin3 <- NULL
for(i in seq(0, 90, 10)) polin1 <- c(polin1, pol1+i)
for(i in seq(0, 90, 10)) polin2 <- c(polin2, pol2+i)
for(i in seq(0, 90, 10)) polin3 <- c(polin3, pol3+i)

parallelograms <- data.frame(
  x = rep(c(ymin, ymin+10, ymax, ymax-10), 10),
  y = polin1,
  id = rep(1:10, each = 4)
)

triangles <-
  bind_rows(data.frame(x = rep(c(ymin, ymin+10, ymin), 10),
                       y = polin2,
                       id = rep(1:10, each = 3)),
            data.frame(x = rep(c(ymax-10, ymax, ymax), 10),
                       y = polin3,
                       id = rep(11:20, each = 3)))

ggplot()+
  # Lexis grid: periods, ages and cohorts
  geom_vline(xintercept = seq(pmin, pmax, 10),
             linewidth = 0.2, linetype = "dashed", alpha = 0.8, color = "grey30")+
  geom_hline(yintercept = seq(amin, amax, 10),
             linewidth = 0.2, linetype = "dashed", alpha = 0.8, color = "grey30")+
  geom_abline(intercept = seq(-pmax, -(pmin-amax), 10), slope = 1,
              linetype = "dashed", color = "grey30", linewidth = .2, alpha = 0.8)+
  geom_polygon(data = parallelograms, aes(x = x, y = y, group = id),
               fill = "blue", alpha = 0.2, color = "black", linewidth = .3)+
  geom_polygon(data = triangles, aes(x = x, y = y, group = id),
               fill = "red", alpha = 0.2, color = "black", linewidth = .3)+
  # white margin with cohort labels
  geom_rect(aes(xmin = 2030.5, xmax = 2040, ymin = -.5, ymax = 120), fill = "white")+
  geom_text(aes(x = 2032, y = seq(22, 102, 10), label = rev(seq(1930, 2010, 10))),
            angle = 45, hjust = 0, size = 3)+
  geom_text(aes(x = 2031, y = 15, label = "Cohort"),
            angle = 45, hjust = 0, size = 3.5)+
  scale_x_continuous(breaks = seq(1910, 2020, 10))+
  scale_y_continuous(breaks = seq(0, 100, 10))+
  coord_fixed(expand = 0, xlim = c(1930, 2040), ylim = c(0, 110))+
  labs(y = "Age", x = "Period")+
  theme_minimal()

ggsave("figures/manuscript/fig02_lexis_diagram.png",
       w = 5, h = 5, dpi = 600)
