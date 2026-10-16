# ==============================================================================
# 03. Figure 1. Life expectancy at birth (1908-2019) in Spain and the four
#     reference countries (panel A) and the gap between Spain and their
#     unweighted mean (panel B). Shaded bands: Civil War (1936-1939) and the
#     1940-1960 catch-up period; annotated points: gap in 1939, 1940 and 1960.
# Input : data_inter/e0_hmd.rds (01_hmd_life_tables.R)
# Output: figures/manuscript/fig01_life_expectancy_gap.png
# ==============================================================================
rm(list = ls()); gc()
library(tidyverse)
library(patchwork)

e0 <- read_rds("data_inter/e0_hmd.rds")

cts_ref <- c("NLD", "NOR", "SWE", "CHE")
ymin <- 1908 # first year with Spanish data in the HMD
ymax <- 2019

names_cts <- c(ESP = "Spain", NLD = "Netherlands", NOR = "Norway",
               SWE = "Sweden", CHE = "Switzerland",
               REF = "Reference mean (unweighted)")
cols <- c(ESP = "#b2182b", NLD = "#e69f00", NOR = "#56b4e9",
          SWE = "#009e73", CHE = "#cc79a7", REF = "black")

le <-
  e0 %>%
  filter(Code %in% c("ESP", cts_ref),
         Year %in% ymin:ymax) %>%
  select(year = Year, code = Code, ex = Total)

ref <-
  le %>%
  filter(code != "ESP") %>%
  reframe(ex = mean(ex), .by = year) %>%
  mutate(code = "REF")

gap <-
  le %>%
  filter(code == "ESP") %>%
  select(year, esp = ex) %>%
  inner_join(ref %>% select(year, ref = ex), by = "year") %>%
  mutate(gap = esp - ref)

# values quoted in the manuscript
gap %>% filter(year %in% c(1908, 1935, 1939, 1940, 1960, 2019))
gap %>% slice_min(gap)

# shaded periods ---------------------------------------------------------------
bands <- tibble(xmin = c(1936, 1940), xmax = c(1939, 1960),
                lab = c("Civil\nWar", "1940-1960"),
                fill = c("grey55", "grey80"))

# panel A: life expectancy by country ------------------------------------------
dtA <-
  bind_rows(le, ref) %>%
  mutate(code = factor(code, levels = names(names_cts)))

pA <-
  ggplot()+
  geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = bands$fill, alpha = .4)+
  geom_text(data = bands, aes(x = (xmin + xmax) / 2, y = 89, label = lab),
            size = 3, colour = "grey30")+
  geom_line(data = dtA, aes(year, ex, colour = code, linetype = code, linewidth = code))+
  scale_colour_manual(values = cols, labels = names_cts, name = NULL)+
  scale_linetype_manual(values = c(ESP = "solid", NLD = "solid", NOR = "solid",
                                   SWE = "solid", CHE = "solid", REF = "dashed"),
                        labels = names_cts, name = NULL)+
  scale_linewidth_manual(values = c(ESP = 1.2, NLD = .5, NOR = .5,
                                    SWE = .5, CHE = .5, REF = .7),
                         labels = names_cts, name = NULL)+
  scale_x_continuous(breaks = seq(1910, 2010, 10))+
  scale_y_continuous(breaks = seq(30, 90, 10))+
  coord_cartesian(xlim = c(ymin, ymax), ylim = c(28, 92), expand = FALSE)+
  labs(y = "Life expectancy at birth (years)", x = NULL)+
  theme_bw()+
  theme(legend.position = "inside",
        legend.position.inside = c(.8, .3),
        legend.background = element_rect(colour = "grey70"),
        legend.key.width = unit(1.2, "cm"),
        panel.grid.minor = element_blank())

# panel B: gap -----------------------------------------------------------------
pts <-
  gap %>%
  filter(year %in% c(1939, 1940, 1960)) %>%
  mutate(lab = sprintf("%d: %.1f", year, gap),
         xlab = c(1926, 1958, 1968),
         ylab = c(-22, -21, -8))

pB <-
  ggplot()+
  geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = bands$fill, alpha = .4)+
  geom_hline(yintercept = 0, colour = "grey40", linewidth = .4)+
  geom_line(data = gap, aes(year, gap), colour = cols[["ESP"]], linewidth = 1.2)+
  geom_segment(data = pts, aes(x = year, y = gap, xend = xlab, yend = ylab),
               colour = "grey50", linewidth = .3)+
  geom_point(data = pts, aes(year, gap), colour = cols[["ESP"]], size = 2.2)+
  geom_label(data = pts, aes(xlab, ylab, label = lab), size = 3,
             linewidth = 0, fill = "white", colour = "grey20")+
  scale_x_continuous(breaks = seq(1910, 2010, 10))+
  coord_cartesian(xlim = c(ymin, ymax), ylim = c(-24, 3), expand = FALSE)+
  labs(y = "Gap (years)\nSpain - reference mean", x = "Year")+
  theme_bw()+
  theme(panel.grid.minor = element_blank())

(pA / pB) +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

ggsave("figures/manuscript/fig01_life_expectancy_gap.png",
       w = 8, h = 7.5, dpi = 300)
