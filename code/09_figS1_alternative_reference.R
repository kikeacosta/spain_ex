# ==============================================================================
# 09. Supplementary Material S1. Robustness check: Figure 5 replicated with an
#     alternative set of reference countries (Denmark, the United Kingdom,
#     Canada and Belgium), and Table S1.1: gap in life expectancy at birth
#     between each of them and Spain (REF - ESP, both sexes) in 1940, 1950
#     and 1960, and reduction of the 1940 gap by 1960.
#
# The United Kingdom series (HMD code GBR_NP) ends in 2022, so its
# decomposition covers 1940-2022; the other three cover 1940-2023.
#
# Input : data_inter/ltables_hmd.rds, data_inter/e0_hmd.rds (01_hmd_life_tables.R)
# Output: data_inter/tableS1_1_altref_gap_1940_1960.csv
#         data_inter/cohort_conts_gap_changes_altref.RData
#         figures/supplement/figS1_cohort_contributions_gap_altref.png
# Note  : the decomposition loop takes a few minutes (see 07).
# ==============================================================================
rm(list = ls()); gc()
source("code/00_functions.R")

cts_alt <- c("DNK", "GBR_NP", "CAN", "BEL")
names_alt <- c(BEL = "Belgium", CAN = "Canada", DNK = "Denmark",
               GBR_NP = "United Kingdom")
ymin <- 1940
ymax <- 2023

hmd_lt <- read_rds("data_inter/ltables_hmd.rds")
missing <- setdiff(cts_alt, unique(hmd_lt$code))
if (length(missing) > 0) stop("Not in data_inter/ltables_hmd.rds: ", paste(missing, collapse = ", "),
                              ". Add them to `cds` in 01_hmd_life_tables.R and re-download.")
hmd_lt <- hmd_lt %>% filter(code %in% c("ESP", cts_alt))

# last year of each comparison: ymax, unless a series ends earlier (GBR_NP: 2022)
ymax_esp <- max(hmd_lt$year[hmd_lt$code == "ESP"])
ymax_alt <- sapply(cts_alt, function(ct) min(ymax, ymax_esp, max(hmd_lt$year[hmd_lt$code == ct])))
print(ymax_alt)

# Table S1.1: gap at birth (REF - ESP, both sexes) in 1940, 1950 and 1960 ------
e0 <- read_rds("data_inter/e0_hmd.rds")

tab <-
  e0 %>%
  filter(Code %in% c("ESP", cts_alt), Year %in% c(1940, 1950, 1960)) %>%
  select(year = Year, code = Code, ex = Total) %>%
  pivot_wider(names_from = code, values_from = ex) %>%
  pivot_longer(all_of(cts_alt), names_to = "ref", values_to = "ex_ref") %>%
  mutate(gap = ex_ref - ESP) %>%
  select(ref, year, ex_esp = ESP, ex_ref, gap) %>%
  pivot_wider(names_from = year, values_from = c(ex_esp, ex_ref, gap)) %>%
  # relative reduction of the 1940 gap by 1960 (%)
  mutate(gap_red_1940_60 = 100 * (1 - gap_1960 / gap_1940),
         country = names_alt[ref]) %>%
  arrange(country) %>%
  select(country, ref, everything()) %>%
  mutate(across(where(is.numeric), ~ round(.x, 1)))

print(as.data.frame(tab))

write_csv(tab, "data_inter/tableS1_1_altref_gap_1940_1960.csv")

# decomposition (Figure S1.1) --------------------------------------------------
res <- decompose_gap_by_age_groups(c1 = "ESP", cts_ref = cts_alt,
                                   ymin = ymin, ymax = ymax_alt)

dt_gaps2_alt <- res$cohs
dt_gaps_summ_alt <- res$summ

save(dt_gaps2_alt, dt_gaps_summ_alt,
     file = "data_inter/cohort_conts_gap_changes_altref.RData")

plot_scc(dt_gaps2_alt, labels = names_alt)

dir.create("figures/supplement", showWarnings = FALSE)
ggsave("figures/supplement/figS1_cohort_contributions_gap_altref.png",
       w = 6, h = 6, dpi = 600)
