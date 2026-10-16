# ==============================================================================
# 07. Figure 5. Scaled cohort-specific contributions (SCC) to the change in the
#     partial life expectancy gap between Spain and each reference country
#     (CHE, NLD, NOR, SWE), 1940-2019, by sex and 10-year age group.
#
# For each age group the gap in partial life expectancy (reference minus
# Spain) is decomposed by single age every year with the symmetric Arriaga
# method; the change in each age contribution from year y to y + 1 is assigned
# to the cohort aged x in year y and cumulated along the cohort's Lexis
# diagonal (est_coh_cont_diff_countries in 00_functions.R). Within each
# sex x age group x country, contributions are divided by the largest absolute
# contribution (SCC in [-1, 1]). Dark marks at the top of each panel: cohorts
# that crossed the age group at some point during 1940-1960.
#
# Input : data_inter/ltables_hmd.rds (01_hmd_life_tables.R)
# Output: data_inter/cohort_conts_gap_changes.RData (decomposition results)
#         figures/manuscript/fig05_cohort_contributions_gap.png
# Note  : the decomposition loop takes a few minutes (11 age groups x 4
#         countries x 80 years x 3 sexes, in both directions).
# ==============================================================================
rm(list = ls()); gc()
source("code/00_functions.R")

cts_ref <- c("SWE", "CHE", "NOR", "NLD")
ymin <- 1940
ymax <- 2019

hmd_lt <-
  read_rds("data_inter/ltables_hmd.rds") %>%
  filter(code %in% c("ESP", cts_ref))

res <- decompose_gap_by_age_groups(c1 = "ESP", cts_ref = cts_ref,
                                   ymin = ymin, ymax = ymax)

dt_gaps2 <- res$cohs
dt_gaps_summ <- res$summ

save(dt_gaps2, dt_gaps_summ, file = "data_inter/cohort_conts_gap_changes.RData")

# figure -----------------------------------------------------------------------
# load("data_inter/cohort_conts_gap_changes.RData")

plot_scc(dt_gaps2)

ggsave("figures/manuscript/fig05_cohort_contributions_gap.png",
       w = 6, h = 6, dpi = 600)
