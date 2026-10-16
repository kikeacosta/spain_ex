# ==============================================================================
# 05. Figure 3, Supplementary Figure S2 and Supplementary Table S2.1.
#     Detrended cohort effects (rate ratios) on all-cause mortality from
#     Age-drift-Period-Cohort models fitted to single-year death counts and
#     exposures (Epi::apc.fit, factor model, "AdPC" parametrization, drift
#     extracted with the unweighted inner product), 1940-2019, by sex, for
#     ages 0-49, 50-99 and 0-99, in Spain and the four reference countries.
#
#     Figure 3  : Spain, all ages (0-99), with 95% confidence bands.
#     Figure S2 : the five countries, rows = age groups (0-49, 50-99, 0-99).
#     Table S2.1: location and depth of the local minimum of the 1940s cohorts
#                 and of the local maximum of the 1960s cohorts, by country,
#                 sex and age group. Because the drift extraction only tilts
#                 the detrended curve linearly, the depth is measured
#                 relative to the straight line joining the neighbouring
#                 turning points of the curve (prominence, in log RR), which
#                 is invariant to that tilt. The deviation from the mean of the
#                 +/-10 neighbouring cohorts is also reported.
#
# Input : data_inter/deaths_exposures_hmd.rds (02_hmd_deaths_exposures.R)
# Output: figures/manuscript/fig03_cohort_effects_apc.png
#         figures/supplement/figS2_cohort_effects_apc_countries.png
#         data_inter/apc_cohort_effects.rds (all estimated cohort effects)
#         data_inter/tableS2_1_cohort_extrema.csv
# Note  : 30 single-year factor models; takes 15-20 minutes.
# ==============================================================================
rm(list = ls()); gc()
library(tidyverse)
library(Epi)

hmd <- read_rds("data_inter/deaths_exposures_hmd.rds")

est_apc <- function(cd = "ESP", sx = "female",
                    amin = 0, amax = 99,
                    pmin = 1940, pmax = 2019) {
  dt <-
    hmd %>%
    filter(code == cd, sex == sx, age %in% amin:amax, year %in% pmin:pmax) %>%
    rename(D = dts, Y = exposure, A = age, P = year) %>%
    mutate(C = P - A)

  # Poisson APC model with age, period and cohort as factors; drift extracted
  # with the unweighted (unit) inner product; period and cohort effects
  # detrended (zero mean and zero slope)
  apc_c <-
    apc.fit(dt, model = "factor", dr.extr = "1", parm = "AdPC", scale = 10^5,
            print.AOV = FALSE)

  list(age = apc_c$Age %>% as_tibble() %>% mutate(code = cd, sex = sx),
       per = apc_c$Per %>% as_tibble() %>% mutate(code = cd, sex = sx),
       coh = apc_c$Coh %>% as_tibble() %>% mutate(code = cd, sex = sx),
       drift = apc_c$Drift[1, 1])
}

cts <- c("ESP", "SWE", "NOR", "NLD", "CHE")
names_cts <- c(ESP = "Spain", SWE = "Sweden", NOR = "Norway",
               NLD = "Netherlands", CHE = "Switzerland")
pmin <- 1940
pmax <- 2019
age_groups <- list(c(0, 49), c(50, 99), c(0, 99))

# models -----------------------------------------------------------------------
# set refit <- FALSE to reuse the saved effects (data_inter/apc_cohort_effects.rds)
# and only redo the figures and the table
refit <- TRUE

if (refit) {
  cohs <- tibble()
  for (cd in cts) {
    for (ag in age_groups) {
      for (sx in c("female", "male")) {
        cat(paste0(cd, " ", sx, " ages ", ag[1], "-", ag[2], "\n"))
        res <- est_apc(cd = cd, sx = sx, amin = ag[1], amax = ag[2],
                       pmin = pmin, pmax = pmax)
        cohs <-
          bind_rows(cohs,
                    res$coh %>% mutate(ages = paste0(ag[1], "-", ag[2]),
                                       drift = res$drift))
      }
    }
  }

  cohs <-
    cohs %>%
    mutate(ages = factor(ages, levels = c("0-49", "50-99", "0-99")),
           sex = str_to_title(sex),
           country = factor(names_cts[code], levels = names_cts),
           lrr = log(`C-RR`))

  write_rds(cohs, "data_inter/apc_cohort_effects.rds")
} else {
  cohs <- read_rds("data_inter/apc_cohort_effects.rds")
}

# Figure 3: Spain, all ages ----------------------------------------------------
cols_sex <- c("#390099", "#5c8001")

cohs %>%
  filter(code == "ESP", ages == "0-99") %>%
  ggplot() +
  geom_ribbon(aes(Coh, ymin = `2.5%`, ymax = `97.5%`, fill = sex), alpha = 0.3) +
  geom_line(aes(Coh, `C-RR`, col = sex)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_log10(breaks = c(0.2, 0.5, 0.7, 1, 1.5, 2, 4, 5)) +
  scale_x_continuous(breaks = seq(1800, 2020, 20)) +
  scale_color_manual(values = cols_sex) +
  scale_fill_manual(values = cols_sex) +
  labs(y = "Relative Risk", x = "Cohort", color = "Sex", fill = "Sex") +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

ggsave("figures/manuscript/fig03_cohort_effects_apc.png",
       w = 8, h = 4.5)

# Figure S2: five countries, three age groups ----------------------------------
cols_cts <- c(Netherlands = "#e69f00", Norway = "#56b4e9", Sweden = "#009e73",
              Switzerland = "#cc79a7", Spain = "#b2182b")
# cohorts discussed in the text: 1940s advantage and early-1960s disadvantage
bands <- tibble(xmin = c(1940, 1957), xmax = c(1950, 1965))

cohs %>%
  mutate(country = fct_relevel(country, "Spain", after = Inf)) %>%
  ggplot() +
  facet_grid(ages ~ sex) +
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
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

dir.create("figures/supplement", showWarnings = FALSE)
ggsave("figures/supplement/figS2_cohort_effects_apc_countries.png",
       w = 8, h = 8, dpi = 300)

# Table S2.1: local extrema of the detrended cohort curves ---------------------
# local extrema within a window of +/- k cohorts
find_extrema <- function(lrr, k = 5){
  n <- length(lrr)
  idx <- seq_len(n)
  is_max <- sapply(idx, function(i) { j <- max(1, i - k):min(n, i + k); lrr[i] == max(lrr[j]) })
  is_min <- sapply(idx, function(i) { j <- max(1, i - k):min(n, i + k); lrr[i] == min(lrr[j]) })
  list(is_max = is_max, is_min = is_min)
}

# depth of the most extreme local minimum (type = "min") or maximum ("max")
# with cohort in `window`, relative to the line joining the nearest turning
# points of opposite type on each side (prominence)
extremum_depth <- function(coh, lrr, window, type = c("min", "max"), k = 5){
  type <- match.arg(type)
  ex <- find_extrema(lrr, k)
  n <- length(lrr)
  cand <- which((if (type == "min") ex$is_min else ex$is_max) & coh %in% window)
  if (length(cand) == 0)
    return(tibble(cohort = NA, lrr = NA, prominence = NA, dev_nb10 = NA,
                  turn_left = NA, turn_right = NA))
  i <- cand[if (type == "min") which.min(lrr[cand]) else which.max(lrr[cand])]
  opp <- which(if (type == "min") ex$is_max else ex$is_min)
  left  <- if (any(opp < i)) max(opp[opp < i]) else 1
  right <- if (any(opp > i)) min(opp[opp > i]) else n
  base <- lrr[left] + (lrr[right] - lrr[left]) * (i - left) / (right - left)
  nb <- setdiff(max(1, i - 10):min(n, i + 10), i)
  # (computed before tibble() so that the new `lrr` column does not shadow the vector)
  lrr_i <- lrr[i]
  prominence <- lrr_i - base
  dev_nb10 <- lrr_i - mean(lrr[nb])
  tibble(cohort = coh[i], lrr = lrr_i, prominence = prominence,
         dev_nb10 = dev_nb10,
         turn_left = coh[left], turn_right = coh[right])
}

features <- tribble(
  ~feature,            ~type, ~window,
  "1940s minimum",     "min", 1935:1955,
  "1960s maximum",     "max", 1955:1970
)

tab <-
  cohs %>%
  arrange(country, sex, ages, Coh) %>%
  group_by(country, sex, ages) %>%
  group_modify(function(d, g) {
    map_dfr(seq_len(nrow(features)), function(f) {
      extremum_depth(d$Coh, d$lrr, features$window[[f]], features$type[f]) %>%
        mutate(feature = features$feature[f], .before = 1)
    })
  }) %>%
  ungroup() %>%
  mutate(across(c(lrr, prominence, dev_nb10), ~ round(.x, 3))) %>%
  arrange(feature, ages, sex, country)

print(as.data.frame(tab))

write_csv(tab, "data_inter/tableS2_1_cohort_extrema.csv")
