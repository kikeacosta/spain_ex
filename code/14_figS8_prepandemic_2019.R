# ==============================================================================
# 14. Supplementary Material S8. Results for the pre-pandemic period 1940-2019
#     (the period of the original submission) compared with 1940-2023.
#
#     Figure S8.1: detrended cohort effects for Spain estimated with data for
#                  1940-2019 and 1940-2023, by sex and age group (0-49, 50-99,
#                  0-99), same specification as Figure 3.
#     Figure S8.2: Figure 5 with data ending in 2019.
#     Table S8.1 : main indicators in both periods: rate ratio of the 1961
#                  cohort (ages 0-99), prominence of the local extremes of
#                  Table S2.1 (Spain and range of the reference countries),
#                  number of cohort contributions tested, % exceptional,
#                  replicated single peaks, opposite-signed pairs found and
#                  replicated, and the replicated peaks and pairs with the
#                  alternative reference set.
#
# The 1940-2023 results are read from the outputs of scripts 05, 07 and 09;
# the 1940-2019 results are re-estimated here with the same functions
# (00_functions.R): 30 APC factor models (15-20 minutes) and two
# decompositions (a few minutes each).
#
# Input : data_inter/deaths_exposures_hmd.rds, data_inter/ltables_hmd.rds
#         data_inter/apc_cohort_effects.rds (05_fig3_cohort_effects_apc.R)
#         data_inter/cohort_conts_gap_changes.RData (07_fig5_cohort_contributions_gap.R)
#         data_inter/cohort_conts_gap_changes_altref.RData (09_figS1_alternative_reference.R)
# Output: data_inter/apc_cohort_effects_2019.rds
#         data_inter/tableS2_1_cohort_extrema_2019.csv
#         data_inter/cohort_conts_gap_changes_2019.RData
#         data_inter/cohort_conts_gap_changes_altref_2019.RData
#         data_inter/tableS8_1_comparison_2019_2023.csv
#         figures/supplement/figS8_1_cohort_effects_apc_2019_vs_2023.png
#         figures/supplement/figS8_2_cohort_contributions_gap_2019.png
# ==============================================================================
rm(list = ls()); gc()
source("code/00_functions.R")
library(Epi)

pmin <- 1940
pmax_old <- 2019
pmax_new <- 2023
per_old <- paste0(pmin, "-", pmax_old)
per_new <- paste0(pmin, "-", pmax_new)

cts <- c("ESP", "SWE", "NOR", "NLD", "CHE")
names_cts <- c(ESP = "Spain", SWE = "Sweden", NOR = "Norway",
               NLD = "Netherlands", CHE = "Switzerland")
cts_ref <- c("SWE", "CHE", "NOR", "NLD")
cts_alt <- c("DNK", "GBR_NP", "CAN", "BEL")

# criterion parameters (as in 11_exceptional_cohorts.R)
k_each_side <- 5; z_thr <- 3; min_ags <- 2; min_cts <- 3; L_max <- 15

# 1. detrended cohort effects with data to 2019 --------------------------------
# set refit <- FALSE to reuse the saved effects (data_inter/apc_cohort_effects_2019.rds)
refit <- TRUE
hmd <- read_rds("data_inter/deaths_exposures_hmd.rds")

if (refit) {
  cohs_old <- fit_cohort_effects(hmd, cts, names_cts, pmin = pmin, pmax = pmax_old)
  write_rds(cohs_old, "data_inter/apc_cohort_effects_2019.rds")
} else {
  cohs_old <- read_rds("data_inter/apc_cohort_effects_2019.rds")
}
rm(hmd)

cohs_new <-
  read_rds("data_inter/apc_cohort_effects.rds") %>%
  mutate(period = per_new)

tab_old <- cohort_extrema_table(cohs_old)
write_csv(tab_old, "data_inter/tableS2_1_cohort_extrema_2019.csv")
tab_new <- cohort_extrema_table(cohs_new)  # = data_inter/tableS2_1_cohort_extrema.csv

# Figure S8.1: Spain, both periods overlaid
cols_per <- setNames(c("grey55", "#b2182b"), c(per_old, per_new))

bind_rows(cohs_old, cohs_new) %>%
  filter(code == "ESP") %>%
  mutate(ages = factor(ages, levels = c("0-49", "50-99", "0-99")),
         period = factor(period, levels = names(cols_per))) %>%
  ggplot() +
  facet_grid(ages ~ sex) +
  geom_line(aes(Coh, `C-RR`, col = period, linewidth = period)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_log10(breaks = c(0.2, 0.5, 0.7, 1, 1.5, 2, 4, 5)) +
  scale_x_continuous(breaks = seq(1850, 2010, 20)) +
  scale_color_manual(values = cols_per) +
  scale_linewidth_manual(values = setNames(c(1, .5), names(cols_per)), guide = "none") +
  labs(y = "Relative Risk", x = "Cohort", color = "Data") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

dir.create("figures/supplement", showWarnings = FALSE)
ggsave("figures/supplement/figS8_1_cohort_effects_apc_2019_vs_2023.png",
       w = 8, h = 8, dpi = 300)

# 2. decompositions with data to 2019 (Figure S8.2) ----------------------------
hmd_lt_all <- read_rds("data_inter/ltables_hmd.rds")

# main reference set (est_age_cont() reads the life tables from `hmd_lt`)
hmd_lt <- hmd_lt_all %>% filter(code %in% c("ESP", cts_ref))
res_old <- decompose_gap_by_age_groups(c1 = "ESP", cts_ref = cts_ref,
                                       ymin = pmin, ymax = pmax_old)
dt_gaps2_2019 <- res_old$cohs
dt_gaps_summ_2019 <- res_old$summ
save(dt_gaps2_2019, dt_gaps_summ_2019,
     file = "data_inter/cohort_conts_gap_changes_2019.RData")

plot_scc(dt_gaps2_2019)
ggsave("figures/supplement/figS8_2_cohort_contributions_gap_2019.png",
       w = 6, h = 6, dpi = 600)

# alternative reference set
hmd_lt <- hmd_lt_all %>% filter(code %in% c("ESP", cts_alt))
res_alt_old <- decompose_gap_by_age_groups(c1 = "ESP", cts_ref = cts_alt,
                                           ymin = pmin, ymax = pmax_old)
dt_gaps2_alt_2019 <- res_alt_old$cohs
dt_gaps_summ_alt_2019 <- res_alt_old$summ
save(dt_gaps2_alt_2019, dt_gaps_summ_alt_2019,
     file = "data_inter/cohort_conts_gap_changes_altref_2019.RData")

# 1940-2023 results of scripts 07 and 09
load("data_inter/cohort_conts_gap_changes.RData")         # dt_gaps2
load("data_inter/cohort_conts_gap_changes_altref.RData")  # dt_gaps2_alt

# 3. exceptionality criterion in both periods ----------------------------------
crit <- function(d) exceptionality_criterion(d, k = k_each_side, z = z_thr,
                                             min_cts = min_cts, min_ags = min_ags,
                                             L_max = L_max)
ex_old     <- crit(dt_gaps2_2019)
ex_new     <- crit(dt_gaps2)
ex_alt_old <- crit(dt_gaps2_alt_2019)
ex_alt_new <- crit(dt_gaps2_alt)

# 4. Table S8.1 ----------------------------------------------------------------
# prominence of an extreme in Spain and its range over the reference countries
prom <- function(tab, feat, sx, ag){
  d <- tab %>% filter(feature == feat, sex == sx, ages == ag) %>%
    mutate(country = as.character(country))
  esp <- d %>% filter(country == "Spain") %>% pull(prominence)
  ref <- d %>% filter(country != "Spain", !is.na(prominence))
  rng <- if (nrow(ref) == 0) "no local extreme in the reference countries" else
    if (nrow(ref) == 1) sprintf("%+.2f, %s only", ref$prominence, ref$country) else
      sprintf("%+.2f to %+.2f", min(ref$prominence), max(ref$prominence))
  sprintf("%+.2f (%s)", esp, rng)
}

# rate ratio of the 1961 cohort, all ages, men / women
rr1961 <- function(cohs){
  d <- cohs %>% filter(code == "ESP", ages == "0-99", Coh == 1961)
  sprintf("%.2f / %.2f", d$`C-RR`[d$sex == "Male"], d$`C-RR`[d$sex == "Female"])
}

summarise_period <- function(cohs, tab, ex, ex_alt){
  c(rr1961(cohs),
    prom(tab, "1940s minimum", "Female", "50-99"),
    prom(tab, "1960s maximum", "Female", "0-49"),
    prom(tab, "1960s maximum", "Male", "0-49"),
    prom(tab, "1960s maximum", "Male", "0-99"),
    format(ex$n_tested, big.mark = ","),
    sprintf("%.1f", ex$pct_exceptional),
    format_peaks(ex$peaks),
    sprintf("%d / %d", nrow(ex$pairs), nrow(ex$blocks)),
    sprintf("%s / %s", format_peaks(ex_alt$peaks), format_blocks(ex_alt$blocks)))
}

tab_s8 <- tibble(
  analysis = c("d-APC (Figure 3, ages 0-99)",
               rep("d-APC (Table S2.1)", 4),
               rep("Decomposition (Figure 5)", 4),
               "Alternative reference set (S1)"),
  indicator = c("RR of the 1961 cohort, men / women",
                "1940s minimum, women 50-99: Spain (reference range)",
                "1960s maximum, women 0-49: Spain (reference range)",
                "1960s maximum, men 0-49: Spain (reference range)",
                "1960s maximum, men 0-99: Spain (reference range)",
                "Cohort contributions tested (ages 0-99)",
                "Exceptional contributions (%)",
                "Single peaks meeting the replication requirement",
                "Opposite-signed pairs found / replicated",
                "Single peaks / pairs meeting the replication requirement"),
  period_1940_2019 = summarise_period(cohs_old, tab_old, ex_old, ex_alt_old),
  period_1940_2023 = summarise_period(cohs_new, tab_new, ex_new, ex_alt_new))

print(as.data.frame(tab_s8))

write_csv(tab_s8, "data_inter/tableS8_1_comparison_2019_2023.csv")
