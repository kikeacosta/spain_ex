# ==============================================================================
# 13. Supplementary Material S7. Spanish life expectancy at birth against its
#     pre-war trend: a linear trend fitted to 1920-1935, by sex (both sexes,
#     women, men), is extrapolated to 1960 and compared with the observed
#     values.
#
#     Figure S7.1: observed life expectancy (solid) and the extrapolated trend
#                  (dashed), 1920-1960, with the shortfall shaded and the
#                  Civil War (1936-1939) marked.
#     Table S7.1 : observed value, trend and deviation (years) in 1935, 1939,
#                  1941, 1945, 1950, 1955 and 1960.
#
# Input : data_inter/e0_hmd.rds (01_hmd_life_tables.R)
# Output: data_inter/tableS7_1_e0_vs_prewar_trend.csv
#         figures/supplement/figS7_1_e0_prewar_trend.png
# ==============================================================================
rm(list = ls()); gc()
library(tidyverse)

e0 <- read_rds("data_inter/e0_hmd.rds")

fit_years <- 1920:1935   # years the trend is fitted to
ymin <- 1920
ymax <- 1960             # last year of the extrapolation
table_years <- c(1935, 1939, 1941, 1945, 1950, 1955, 1960)

le <-
  e0 %>%
  filter(Code == "ESP", Year %in% ymin:ymax) %>%
  pivot_longer(c(Total, Female, Male), names_to = "sex", values_to = "ex") %>%
  mutate(sex = factor(sex, levels = c("Total", "Female", "Male"),
                      labels = c("Both sexes", "Women", "Men"))) %>%
  select(sex, year = Year, ex)

# linear trend fitted to 1920-1935 and extrapolated to 1960, by sex -------------
trend <-
  le %>%
  group_by(sex) %>%
  group_modify(function(d, g){
    m <- lm(ex ~ year, data = d %>% filter(year %in% fit_years))
    d %>% mutate(trend = predict(m, newdata = d),
                 slope = coef(m)[["year"]])
  }) %>%
  ungroup() %>%
  mutate(deviation = ex - trend)

# slope of the trend (years of life expectancy gained per calendar year)
print(trend %>% distinct(sex, slope) %>% mutate(slope = round(slope, 3)))

# largest shortfall
print(trend %>% group_by(sex) %>% slice_min(deviation, n = 1) %>% ungroup() %>%
        mutate(across(c(ex, trend, deviation), ~ round(.x, 1))))

# Table S7.1 -------------------------------------------------------------------
tab <-
  trend %>%
  filter(year %in% table_years) %>%
  select(sex, year, observed = ex, trend, deviation) %>%
  mutate(across(c(observed, trend, deviation), ~ round(.x, 1))) %>%
  arrange(sex, year)

print(as.data.frame(tab))

write_csv(tab, "data_inter/tableS7_1_e0_vs_prewar_trend.csv")

# Figure S7.1 ------------------------------------------------------------------
war <- tibble(xmin = 1936, xmax = 1939)

trend %>%
  ggplot() +
  facet_wrap(~ sex) +
  geom_rect(data = war, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey55", alpha = .4) +
  # shortfall: observed below the trend
  geom_ribbon(aes(year, ymin = pmin(ex, trend), ymax = trend),
              fill = "#b2182b", alpha = .25) +
  geom_line(aes(year, trend), linetype = "dashed", linewidth = .5, colour = "grey20") +
  geom_line(aes(year, ex), linewidth = .9, colour = "#b2182b") +
  scale_x_continuous(breaks = seq(1920, 1960, 10)) +
  scale_y_continuous(breaks = seq(30, 80, 10)) +
  coord_cartesian(xlim = c(ymin - 1, ymax + 1), expand = FALSE) +
  labs(x = "Year", y = "Life expectancy at birth (years)") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(.8, "cm"))

dir.create("figures/supplement", showWarnings = FALSE)
ggsave("figures/supplement/figS7_1_e0_prewar_trend.png", w = 9, h = 3.8, dpi = 300)
