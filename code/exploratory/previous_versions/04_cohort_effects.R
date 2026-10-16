rm(list = ls())
source("code/00_functions.R")
library(Epi)
library(patchwork)
library(tidyverse)

hmd <- read_rds("data_inter/deaths_exposures_hmd.rds")
unique(hmd$sex)
unique(hmd$age)
unique(hmd$year)
unique(hmd$code)

# cd = "USA"; sx = "female";
# amin = 0; amax = 90;
# pmin = 1970; pmax = 2019

est_apc <- function(
  cd = "USA",
  sx = "female",
  amin = 0,
  amax = 90,
  pmin = 1970,
  pmax = 2019
) {
  dt <-
    hmd %>%
    filter(code == cd, sex == sx, age %in% amin:amax, year %in% pmin:pmax) %>%
    rename(D = dts, Y = exposure, A = age, P = year) %>%
    # creating new age categories grouped in 5 years
    mutate(C = P - A)
  # aggregating deaths and exposures by age and period intervals

  # fitting the model
  # ~~~~~~~~~~~~~~~~~
  apc_c <-
    apc.fit(dt, model = "factor", dr.extr = "1", parm = "AdPC", scale = 10^5)

  as <-
    apc_c$Age %>%
    as_tibble() %>%
    mutate(code = cd, sex = sx)

  ps <-
    apc_c$Per %>%
    as_tibble() %>%
    mutate(code = cd, sex = sx)

  cs <-
    apc_c$Coh %>%
    as_tibble() %>%
    mutate(code = cd, sex = sx)

  out <- list(as, ps, cs)

  return(out)
}

{
  cd <- "ESP"
  amin <- 0
  amax <- 99
  pmin <- 1940
  pmax <- 2019
  nplot <- "1d"

  apc_f_y <- est_apc(
    cd = cd,
    sx = "female",
    amin = 0,
    amax = 49,
    pmin = pmin,
    pmax = pmax
  )

  apc_m_y <- est_apc(
    cd = cd,
    sx = "male",
    amin = 0,
    amax = 49,
    pmin = pmin,
    pmax = pmax
  )

  cohs_y <-
    bind_rows(
      apc_f_y[3],
      apc_m_y[3],
    ) %>%
    mutate(ages = "0-49")

  apc_f_o <- est_apc(
    cd = cd,
    sx = "female",
    amin = 50,
    amax = 99,
    pmin = pmin,
    pmax = pmax
  )

  apc_m_o <- est_apc(
    cd = cd,
    sx = "male",
    amin = 50,
    amax = 99,
    pmin = pmin,
    pmax = pmax
  )

  cohs_o <-
    bind_rows(
      apc_f_o[3],
      apc_m_o[3],
    ) %>%
    mutate(ages = "50-99")

  apc_f_t <- est_apc(
    cd = cd,
    sx = "female",
    amin = 0,
    amax = 99,
    pmin = pmin,
    pmax = pmax
  )

  apc_m_t <- est_apc(
    cd = cd,
    sx = "male",
    amin = 0,
    amax = 99,
    pmin = pmin,
    pmax = pmax
  )

  cohs_t <-
    bind_rows(
      apc_f_t[3],
      apc_m_t[3],
    ) %>%
    mutate(ages = "0-99")

  unique(cohs$Coh) %>% sort

  cols <- c("#390099", "#5c8001")

  cohs_t %>%
    # filter(Coh %in% 1880:1980) %>%
    ggplot() +
    geom_ribbon(
      aes(Coh, ymin = `2.5%`, ymax = `97.5%`, fill = sex),
      alpha = 0.3
    ) +
    geom_line(aes(Coh, `C-RR`, col = sex)) +
    # geom_vline(xintercept = c(1900, 1960))+
    scale_y_log10(breaks = c(0.2, 0.5, 0.7, 1, 1.5, 2, 4, 5)) +
    scale_x_continuous(breaks = seq(1800, 2020, 10)) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    theme_bw() +
    coord_cartesian(
      # xlim = c(1900, 2020),
      # ylim = c(1/2.5, 2.1)
    ) +
    labs(
      y = "Relative Risk",
      x = "Cohort",
      color = "Sex",
      fill = "Sex",
      title = paste0("Ages ", amin, "-", amax, ", period ", pmin, "-", pmax),
    ) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(legend.position = "bottom", plot.title = element_text(size = 9))

  ggsave(
    paste0(
      "figures/manuscript/fig0",
      nplot,
      "_cohort_effects_",
      cd,
      "_",
      amin,
      "_",
      amax,
      "_",
      pmin,
      "_",
      pmax,
      ".png"
    ),
    w = 8,
    h = 4
  )
}
# ggsave(paste0("figures/apc_effects/cohort_effects_", cd, "_", amin, "_", amax, "_",
#               pmin, "_", pmax, ".png"),
#        w = 8, h = 4)
cohs <-
  bind_rows(cohs_y, cohs_o, cohs_t) %>%
  mutate(
    ages = factor(ages, levels = c("0-49", "50-99", "0-99")),
    sex = str_to_title(sex)
  )

cols <- c("#390099", "#5c8001")

cohs %>%
  # filter(Coh %in% 1880:1980) %>%
  ggplot() +
  facet_grid(ages ~ .) +
  geom_ribbon(
    aes(Coh, ymin = `2.5%`, ymax = `97.5%`, fill = sex),
    alpha = 0.3
  ) +
  geom_line(aes(Coh, `C-RR`, col = sex)) +
  # geom_vline(xintercept = c(1900, 1960))+
  scale_y_log10(breaks = c(0.2, 0.5, 0.7, 1, 1.5, 2, 4, 5)) +
  scale_x_continuous(breaks = seq(1800, 2020, 10)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  coord_cartesian(
    # xlim = c(1900, 2020),
    # ylim = c(1/2.5, 2.1)
  ) +
  labs(
    y = "Relative Risk",
    x = "Cohort",
    color = "Sex",
    fill = "Sex",
    # title = paste0("Ages ", amin, "-", amax, ", period ", pmin, "-", pmax),
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 9),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(
  "figures/manuscript/fig00_cohort_effects_3_age_groups.png",
  w = 8,
  h = 5.5
)
