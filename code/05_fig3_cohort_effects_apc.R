# ==============================================================================
# 05. Figure 3, Supplementary Figures S2.1-S2.5 and Supplementary Table S2.1.
#     Detrended cohort effects (rate ratios) on all-cause mortality from
#     Age-drift-Period-Cohort models fitted to single-year death counts and
#     exposures (Epi::apc.fit, factor model, "AdPC" parametrization, drift
#     extracted with the unweighted inner product), 1940-2023, by sex, for
#     ages 0-49, 50-99 and 0-99, in Spain and the four reference countries
#     (est_apc() and fit_cohort_effects() in 00_functions.R).
#
#     Figure 3   : Spain, all ages (0-99), with 95% confidence bands.
#     Figure S2.1: reference countries, ages 0-99, with 95% confidence bands.
#     Figure S2.2: the five countries overlaid, ages 0-99, by sex.
#     Figure S2.3: Spain, ages 0-49 and 50-99, with 95% confidence bands.
#     Figure S2.4: reference countries, ages 0-49 and 50-99, with 95% bands.
#     Figure S2.5: the five countries overlaid, ages 0-49 and 50-99, by sex.
#     Table S2.1 : location and depth of the local minimum of the 1940s cohorts
#                  and of the local maximum of the 1960s cohorts, by country,
#                  sex and age group. Because the drift extraction only tilts
#                  the detrended curve linearly, the depth is measured
#                  relative to the straight line joining the neighbouring
#                  turning points of the curve (prominence, in log RR), which
#                  is invariant to that tilt (cohort_extrema_table() in
#                  00_functions.R). The deviation from the mean of the +/-10
#                  neighbouring cohorts is also reported.
#
# Input : data_inter/deaths_exposures_hmd.rds (02_hmd_deaths_exposures.R)
# Output: figures/manuscript/fig03_cohort_effects_apc.png
#         figures/supplement/figS2_1_cohort_effects_apc_ref_0_99.png
#         figures/supplement/figS2_2_cohort_effects_apc_countries_0_99.png
#         figures/supplement/figS2_3_cohort_effects_apc_spain_age_groups.png
#         figures/supplement/figS2_4_cohort_effects_apc_ref_age_groups.png
#         figures/supplement/figS2_5_cohort_effects_apc_countries_age_groups.png
#         data_inter/apc_cohort_effects.rds (all estimated cohort effects)
#         data_inter/tableS2_1_cohort_extrema.csv
# Note  : 30 single-year factor models; takes 15-20 minutes.
# ==============================================================================
rm(list = ls()); gc()
source("code/00_functions.R")
library(Epi)

hmd <- read_rds("data_inter/deaths_exposures_hmd.rds")

cts <- c("ESP", "SWE", "NOR", "NLD", "CHE")
names_cts <- c(ESP = "Spain", SWE = "Sweden", NOR = "Norway",
               NLD = "Netherlands", CHE = "Switzerland")
pmin <- 1940
pmax <- 2023
age_groups <- list(c(0, 49), c(50, 99), c(0, 99))

# models -----------------------------------------------------------------------
# set refit <- FALSE to reuse the saved effects (data_inter/apc_cohort_effects.rds)
# and only redo the figures and the table
refit <- TRUE

if (refit) {
  cohs <- fit_cohort_effects(hmd, cts, names_cts, age_groups, pmin, pmax)
  write_rds(cohs, "data_inter/apc_cohort_effects.rds")
} else {
  cohs <- read_rds("data_inter/apc_cohort_effects.rds")
}

# figures ----------------------------------------------------------------------
cols_sex <- c(Female = "#390099", Male = "#5c8001")
cols_cts <- c(Netherlands = "#e69f00", Norway = "#56b4e9", Sweden = "#009e73",
              Switzerland = "#cc79a7", Spain = "#b2182b")
# cohorts discussed in the text: 1940s advantage and early-1960s disadvantage
bands <- tibble(xmin = c(1940, 1957), xmax = c(1950, 1965))

theme_apc <- function(){
  theme_bw() +
    theme(legend.position = "bottom",
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank())
}

# cohort curves with 95% confidence bands, coloured by sex
plot_ci <- function(d){
  d %>%
    ggplot() +
    geom_ribbon(aes(Coh, ymin = `2.5%`, ymax = `97.5%`, fill = sex), alpha = 0.3) +
    geom_line(aes(Coh, `C-RR`, col = sex)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_log10(breaks = c(0.2, 0.5, 0.7, 1, 1.5, 2, 4, 5)) +
    scale_x_continuous(breaks = seq(1800, 2020, 20)) +
    scale_color_manual(values = cols_sex) +
    scale_fill_manual(values = cols_sex) +
    labs(y = "Relative Risk", x = "Cohort", color = "Sex", fill = "Sex") +
    theme_apc()
}

# the five countries overlaid, Spain highlighted
plot_countries <- function(d){
  d %>%
    mutate(country = fct_relevel(country, "Spain", after = Inf)) %>%
    ggplot() +
    geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "grey85", alpha = .5) +
    geom_line(aes(Coh, `C-RR`, col = country, linewidth = country)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_log10(breaks = c(0.2, 0.5, 0.7, 1, 1.5)) +
    scale_x_continuous(breaks = seq(1850, 2010, 20)) +
    scale_color_manual(values = cols_cts) +
    scale_linewidth_manual(values = c(Netherlands = .4, Norway = .4, Sweden = .4,
                                      Switzerland = .4, Spain = .9), guide = "none") +
    labs(y = "Relative Risk", x = "Cohort", color = NULL) +
    theme_apc()
}

dir.create("figures/supplement", showWarnings = FALSE)

# Figure 3: Spain, all ages
plot_ci(cohs %>% filter(code == "ESP", ages == "0-99"))
ggsave("figures/manuscript/fig03_cohort_effects_apc.png", w = 8, h = 4.5)

# Figure S2.1: reference countries, all ages
plot_ci(cohs %>% filter(code != "ESP", ages == "0-99")) +
  facet_wrap(~ country, ncol = 2)
ggsave("figures/supplement/figS2_1_cohort_effects_apc_ref_0_99.png",
       w = 8, h = 6, dpi = 300)

# Figure S2.2: five countries, all ages, by sex
plot_countries(cohs %>% filter(ages == "0-99")) +
  facet_wrap(~ sex)
ggsave("figures/supplement/figS2_2_cohort_effects_apc_countries_0_99.png",
       w = 8, h = 4.5, dpi = 300)

# Figure S2.3: Spain, ages 0-49 (top) and 50-99 (bottom)
plot_ci(cohs %>% filter(code == "ESP", ages != "0-99")) +
  facet_grid(ages ~ .)
ggsave("figures/supplement/figS2_3_cohort_effects_apc_spain_age_groups.png",
       w = 8, h = 7, dpi = 300)

# Figure S2.4: reference countries, ages 0-49 (left) and 50-99 (right)
plot_ci(cohs %>% filter(code != "ESP", ages != "0-99")) +
  facet_grid(country ~ ages)
ggsave("figures/supplement/figS2_4_cohort_effects_apc_ref_age_groups.png",
       w = 8, h = 10, dpi = 300)

# Figure S2.5: five countries, ages 0-49 (top) and 50-99 (bottom), by sex
plot_countries(cohs %>% filter(ages != "0-99")) +
  facet_grid(ages ~ sex)
ggsave("figures/supplement/figS2_5_cohort_effects_apc_countries_age_groups.png",
       w = 8, h = 6.5, dpi = 300)

# Table S2.1: local extrema of the detrended cohort curves ---------------------
tab <- cohort_extrema_table(cohs)

print(as.data.frame(tab))

write_csv(tab, "data_inter/tableS2_1_cohort_extrema.csv")
