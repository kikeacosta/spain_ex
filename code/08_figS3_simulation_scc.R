# ==============================================================================
# 08. Supplementary Figure S3. Simulation of the cohort-indexed decomposition
#     (Step 3 of the manuscript) under three scenarios for the mortality of a
#     hypothetical Spain relative to a reference population:
#       A. period-only, nonlinear catch-up in 1940-1960
#       B. A plus a block of cohorts (1945-1954) with 15% lower mortality
#       C. period-only, linear convergence 1940-2000
#     The decomposition is run with the same functions used for Figure 5
#     (00_functions.R), so the simulation validates the actual pipeline:
#     a period-only catch-up produces a smooth band over all the cohorts that
#     crossed each age interval during the catch-up years, whereas a block of
#     advantaged cohorts produces opposite-signed spikes at its edges.
#
# R port of simulation_SCC.py (code/exploratory/python/, written by the
# coauthor); same baseline, scenarios and layout.
# Output: figures/supplement/figS3_simulation_scc.png
# ==============================================================================
rm(list = ls()); gc()
source("code/00_functions.R")

# simulated mortality surfaces --------------------------------------------------
ages  <- 0:110
years <- 1900:2019

grid <-
  expand_grid(age = ages, year = years) %>%
  mutate(cohort = year - age)

# Siler-type baseline hazard (reference population c. 1900), declining 2% a year
base <- function(x) 0.12 * exp(-1.2 * x) + 0.004 + 0.00008 * exp(0.09 * x)

mx <-
  grid %>%
  mutate(m_ref = base(age) * exp(-0.02 * (year - 1900)),
         # Spain: excess log-mortality flat until 1940, rapid (nonlinear)
         # catch-up 1940-1960, gone by about 1980
         G = 0.9 - 0.85 / (1 + exp(-(year - 1950) / 4)),
         G = ifelse(year > 1980, G * exp(-(year - 1980) / 5), G),
         # homogeneous (linear) convergence 1940-2000
         lin = pmin(pmax(0.9 * (1 - (year - 1940) / 60), 0), 0.9),
         REF   = m_ref,
         ESP_A = m_ref * exp(G),
         ESP_B = m_ref * exp(G) * ifelse(cohort %in% 1945:1954, 0.85, 1),
         ESP_C = m_ref * exp(lin)) %>%
  select(age, year, REF, ESP_A, ESP_B, ESP_C) %>%
  pivot_longer(c(REF, ESP_A, ESP_B, ESP_C), names_to = "code", values_to = "mx")

# life tables in the format expected by est_age_cont(): code, year, sex, age,
# lx, Lx (constant hazard within each single year of age, radix 1)
hmd_lt <-
  mx %>%
  arrange(code, year, age) %>%
  group_by(code, year) %>%
  mutate(qx = 1 - exp(-mx),
         lx = cumprod(c(1, 1 - qx[-n()])),
         Lx = (lx + lead(lx, default = 0)) / 2) %>%
  ungroup() %>%
  mutate(sex = "t")

# sanity check: the symmetric contributions sum to e2 - e1 (ages 0-9, 1945)
chk <- est_age_cont(c1 = "ESP_A", c2 = "REF", y1 = 1945, y2 = 1945,
                    amin = 0, amax = 9)
cat("check:", sum(chk$dt), "vs", chk$ex2[1] - chk$ex1[1], "\n")

# decomposition with the project functions --------------------------------------
scen <- c(ESP_A = "A. Period only\n(nonlinear catch-up 1940-60)",
          ESP_B = "B. A + cohort advantage\n(cohorts 1945-54, -15%)",
          ESP_C = "C. Period only\n(linear convergence 1940-2000)")

sim <- tibble()
for (a in seq(0, 90, 10)) {
  for (sc in names(scen)) {
    cat(paste0("ages ", a, "-", a + 9, " ", sc, "\n"))
    res <- est_coh_cont_diff_countries(c1 = sc, c2 = "REF",
                                       amin = a, amax = a + 9,
                                       ymin = 1940, ymax = 2019,
                                       sexes = "t")
    sim <- bind_rows(sim, res$coh_cont)
  }
}

sim2 <-
  sim %>%
  mutate(amin = as.integer(str_extract(ags, "^[0-9]+")),
         ags = str_replace(ags, "_", "-"),
         ags = fct_rev(fct_reorder(ags, amin)),
         scen = factor(scen[cts %>% str_remove("_REF")], levels = scen)) %>%
  group_by(scen, ags) %>%
  mutate(scc = dt_coh / max(abs(dt_coh))) %>%
  ungroup()

write_rds(sim2, "data_inter/simulation_scc.rds")

# bands: cohorts that crossed each age interval during 1940-1960 (grey) and
# the advantaged cohorts of scenario B (yellow)
bands <-
  sim2 %>%
  distinct(scen, ags, amin) %>%
  mutate(xmin = 1940 - amin - 9, xmax = 1960 - amin)

# figure ------------------------------------------------------------------------
sim2 %>%
  ggplot() +
  facet_grid(ags ~ scen, switch = "y") +
  geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey70", alpha = .4) +
  geom_rect(data = bands %>% filter(str_starts(scen, "B")) %>% mutate(xmin = 1945, xmax = 1955),
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "gold", alpha = .4) +
  geom_col(aes(cohort, scc, fill = scc < 0), width = 1) +
  geom_hline(yintercept = 0, linewidth = .3) +
  scale_fill_manual(values = c(`TRUE` = "#2a6f97", `FALSE` = "#c1121f"), guide = "none") +
  scale_x_continuous(breaks = seq(1860, 2010, 30)) +
  scale_y_continuous(breaks = c(-1, 0, 1), position = "right") +
  coord_cartesian(xlim = c(1850, 2010), ylim = c(-1.1, 1.1), expand = 0) +
  labs(x = "Birth cohort", y = "Scaled cohort contribution (SCC) by age interval",
       caption = paste0("Negative values (blue) = narrowing of the REF-ESP gap; positive (red) = widening. ",
                        "Grey: cohorts that crossed each age interval during 1940-60.\n",
                        "Yellow: cohorts 1945-54 with 15% lower mortality at all ages (scenario B); ",
                        "with data ending in 2019, these cohorts are not observed over a complete interval beyond ages 60-69.")) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 8),
        strip.text.x = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.title.y.right = element_text(angle = 90),
        panel.grid = element_blank(),
        panel.spacing.x = unit(0.6, "cm"),
        plot.caption = element_text(size = 6.5, hjust = .5))

dir.create("figures/supplement", showWarnings = FALSE)
ggsave("figures/supplement/figS3_simulation_scc.png", w = 7.5, h = 10.5, dpi = 300)
