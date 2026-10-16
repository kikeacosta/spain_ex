# ==============================================================================
# 09. Supplementary Figure S1. Robustness check: Figure 5 replicated with an
#     alternative set of reference countries.
#
# The manuscript mentions Denmark, England and Wales and France. England and
# Wales (HMD code GBRTENW) is not in the default download of
# 01_hmd_life_tables.R, so the United Kingdom (GBR_NP) is used unless
# GBRTENW is available in data_inter/ltables_hmd.rds (add it to `cds` in
# 01_hmd_life_tables.R and re-download to use it).
#
# Input : data_inter/ltables_hmd.rds
# Output: data_inter/cohort_conts_gap_changes_altref.RData
#         figures/supplement/figS1_cohort_contributions_gap_altref.png
# ==============================================================================
rm(list = ls()); gc()
source("code/00_functions.R")

hmd_lt <- read_rds("data_inter/ltables_hmd.rds")

cts_alt <- c("DNK", "GBRTENW", "FRATNP")
if (!"GBRTENW" %in% hmd_lt$code) {
  message("GBRTENW (England and Wales) not in the life tables file; using GBR_NP (United Kingdom) instead.")
  cts_alt <- c("DNK", "GBR_NP", "FRATNP")
}
missing <- setdiff(cts_alt, unique(hmd_lt$code))
if (length(missing) > 0) stop("Not in data_inter/ltables_hmd.rds: ", paste(missing, collapse = ", "),
                              ". Add them to `cds` in 01_hmd_life_tables.R and re-download.")

ymin <- 1940
ymax <- 2019

hmd_lt <- hmd_lt %>% filter(code %in% c("ESP", cts_alt))

res <- decompose_gap_by_age_groups(c1 = "ESP", cts_ref = cts_alt,
                                   ymin = ymin, ymax = ymax)

dt_gaps2_alt <- res$cohs
dt_gaps_summ_alt <- res$summ

save(dt_gaps2_alt, dt_gaps_summ_alt,
     file = "data_inter/cohort_conts_gap_changes_altref.RData")

plot_scc(dt_gaps2_alt, cols = c("#AE2012", "#005F73", "#EE9B00"))

dir.create("figures/supplement", showWarnings = FALSE)
ggsave("figures/supplement/figS1_cohort_contributions_gap_altref.png",
       w = 6, h = 6, dpi = 600)
