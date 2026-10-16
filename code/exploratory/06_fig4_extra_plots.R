# ==============================================================================
# Exploratory plots of the deviations from the Age-drift-Period model
# (rates by age, cohort and period; tile plots). Not in the manuscript.
# Requires data_inter/ap_residuals_ESP.rds produced by
# code/06_fig4_ap_residual_surfaces.R
# ==============================================================================
rm(list = ls()); gc()
library(tidyverse)
curvs <- read_rds("data_inter/ap_residuals_ESP.rds")
sx <- "male"


curvs %>%
  filter(sex == sx,
         age %in% seq(0, 50, 10)) %>%
  select(year, age, starts_with("mx"), -mx_spl) %>%
  gather(starts_with("mx"), key = type, value = mx) %>%
  filter(type %in% c("mx", "mx_ap_spl")) %>% 
  ggplot()+
  facet_wrap(~age, scales = "free")+
  scale_y_log10()+
  geom_line(aes(year, mx, col = type))+
  theme_bw()

curvs %>%
  mutate(cohort = year - age) %>%
  filter(sex == "female",
         # year >= 1970,
         age %in% seq(0, 50, 10)) %>%
  select(year, age, cohort, starts_with("mx"), -mx_spl) %>%
  gather(starts_with("mx"), key = type, value = mx) %>%
  filter(type %in% c("mx", "mx_ap_spl")) %>% 
  ggplot()+
  facet_wrap(~age, scales = "free")+
  scale_y_log10()+
  geom_line(aes(cohort, mx, col = type))+
  theme_bw()

curvs %>%
  mutate(cohort = year - age) %>%
  filter(sex == "female",
         year %in% seq(1940, 1980, 10)) %>%
  select(year, age, cohort, starts_with("mx")) %>%
  gather(starts_with("mx"), key = type, value = mx) %>%
  filter(type %in% c("mx", "mx_ap_spl")) %>% 
  ggplot()+
  facet_wrap(~year, scales = "free")+
  geom_line(aes(age, mx, col = type))+
  scale_y_log10()+
  theme_bw()

# curvs %>%
#   mutate(cohort = year - age) %>%
#   filter(sex == sx,
#          year >= 1970,
#          age %in% seq(20, 80, 10)) %>%
#   select(year, age, cohort, starts_with("mx")) %>%
#   gather(starts_with("mx"), key = type, value = mx) %>%
#   filter(type %in% c("mx_exc2", "mx_exc")) %>%
#   ggplot()+
#   facet_wrap(~age, scales = "free")+
#   geom_line(aes(cohort, mx, col = type))+
#   geom_hline(yintercept = 0, lty = "dashed")+
#   theme_bw()

curvs %>%
  mutate(cohort = year - age) %>%
  filter(sex == sx,
         year >= 1970,
         age %in% seq(20, 80, 10)) %>%
  # select(year, age, cohort, starts_with("mx")) %>%
  # gather(starts_with("mx"), key = type, value = mx) %>%
  # mutate(log_mx = log(mx))+
  # filter(type %in% c("mx_exc2", "mx_exc")) %>%
  ggplot()+
  facet_wrap(~age, scales = "free_x")+
  geom_line(aes(cohort, rr))+
  geom_hline(yintercept = 1, lty = "dashed")+
  theme_bw()

curvs %>%
  mutate(cohort = year - age) %>%
  filter(sex == sx,
         year >= 1960) %>%
  ggplot()+
  geom_tile(aes(year, age, fill = rr2))+
  scale_fill_gradient2(midpoint = 1, low = "blue", mid = "white", high = "red")+
  scale_y_continuous(breaks=seq(20, 90, 10)) +
  # geom_vline(xintercept = seq(1960, 2010, by=10), linetype="dashed", color="grey30", size=.15) +
  # geom_hline(yintercept = seq(20, 90, by=10), linetype="dashed", color="grey30", size=.15) +
  # geom_abline(intercept = seq(-2010, -1860, by=10),
  #             slope = 1, linetype="dashed", color="grey30", size=.15)+
  coord_fixed(expand = 0)+
  theme_bw()


