# ==============================================================================
# 11. Operational criterion for exceptional and disproportionate cohort
#     contributions (Methods, Step 3; Supplementary Material S5), its power
#     check on the simulation of Figure S3 and a sensitivity analysis. The
#     same criterion is applied to the alternative reference set of
#     Supplementary Material S1.
#
# What the criterion looks for. The cohort-indexed decomposition of Figure 5
# assigns the change in the LE gap to cohort positions. A period change spreads
# over all the cohorts that crossed an age interval at the time; a cohort
# effect that is confined to a block of cohorts does not produce a sustained
# contribution but a pair of opposite-signed peaks at the two edges of the
# block (Figure S3, scenario B). The criterion therefore flags DISCONTINUITIES
# between adjacent cohorts, in both directions, and then requires replication.
# It does not detect a cohort advantage that builds up gradually over many
# cohorts: such a pattern is indistinguishable from period change in this
# decomposition and is the object of the d-APC analysis (Figures 3 and S2).
#
# Definitions (parameters below; functions in 00_functions.R).
#   Neighbourhood : the k_each_side cohorts on each side of cohort c, in the
#                   same sex, age interval and comparator (default 5 + 5).
#   Robust score  : z_c = (Cc - median of neighbours) / (1.4826 * MAD of
#                   neighbours). The MAD replaces the standard deviation used
#                   in a first version so that a single large neighbour does
#                   not mask a peak (the standard deviation of ten values is
#                   itself a noisy and easily inflated estimate).
#   Exceptional   : |z_c| > z_thr (default 3, the usual cut-off for robust
#                   z-scores), with sign recorded (narrowing / widening).
#                   The first and last k_each_side cohorts of each age
#                   interval are not tested (one-sided neighbourhoods).
#   Disproportionate (single peaks): the same cohort, sex and sign flagged in
#                   at least min_cts comparators within an age interval, in at
#                   least min_ags age intervals. Period changes fall on
#                   different cohorts in different age intervals, so the
#                   replication across age intervals is the condition that
#                   separates cohort from period signals; the replication
#                   across comparators only removes comparator-specific noise,
#                   because the four gap series share the Spanish data.
#   Block signature (pairs): a narrowing peak at cohort c followed by a
#                   widening peak at c + L, 2 <= L <= L_max (an advantaged
#                   block of L cohorts), or the reverse (a disadvantaged
#                   block). A block is disproportionate when the same
#                   (entry, exit) pair is found in at least min_cts
#                   comparators in at least min_ags age intervals.
#
# Input : data_inter/cohort_conts_gap_changes.RData (07_fig5_cohort_contributions_gap.R)
#         data_inter/cohort_conts_gap_changes_altref.RData (09_figS1_alternative_reference.R)
#         data_inter/simulation_scc.rds (08_figS3_simulation_scc.R)
# Output: data_inter/exceptional_cohorts.csv              all flagged contributions
#         data_inter/disproportionate_cohorts.csv         single peaks meeting the criterion
#         data_inter/block_pairs.csv                      opposite-signed pairs found
#         data_inter/tableS5_1_exceptional_by_age_sex.csv flagged contributions by age interval, sex and sign (Table S5.1)
#         data_inter/exceptional_cohorts_simulation.csv   power check on the simulation (Table S5.2)
#         data_inter/exceptional_cohorts_sensitivity.csv  results for other settings (Table S5.3)
#         data_inter/exceptional_cohorts_altref.csv       flagged contributions, alternative reference set (S1)
# ==============================================================================
rm(list = ls()); gc()
source("code/00_functions.R")

# parameters -------------------------------------------------------------------
k_each_side <- 5    # neighbouring cohorts on each side (10 in total)
z_thr       <- 3    # threshold on the robust z-score
min_ags     <- 2    # age intervals in which the signal must recur
min_cts     <- 3    # comparator countries in which it must recur (out of 4)
L_max       <- 15   # longest block of cohorts considered for the pair signature

# observed data ----------------------------------------------------------------
load("data_inter/cohort_conts_gap_changes.RData")

res <- exceptionality_criterion(dt_gaps2, k = k_each_side, z = z_thr,
                                min_cts = min_cts, min_ags = min_ags, L_max = L_max)

cat("tested contributions:", res$n_tested, " exceptional:", res$n_exceptional,
    sprintf("(%.1f%%)", res$pct_exceptional), "\n")

cat("\nsingle peaks meeting the criterion (same sign in >=", min_cts, "comparators and >=", min_ags, "age intervals):\n")
print(as.data.frame(res$peaks))

cat("\nopposite-signed pairs found in the data:", nrow(res$pairs),
    "; pairs meeting the replication requirement:", nrow(res$blocks), "\n")
if (nrow(res$blocks) > 0) print(as.data.frame(res$blocks))

# Table S5.1: exceptional contributions by age interval, sex and sign
tab_s5_1 <-
  res$exceptional %>%
  count(ags, sex, direction) %>%
  mutate(sex = ifelse(sex == "f", "women", "men")) %>%
  pivot_wider(names_from = c(sex, direction), values_from = n, values_fill = 0) %>%
  select(ags, women_narrowing, women_widening, men_narrowing, men_widening) %>%
  arrange(ags)

cat("\nexceptional contributions by age interval, sex and sign (Table S5.1):\n")
print(as.data.frame(tab_s5_1))

write_csv(res$exceptional, "data_inter/exceptional_cohorts.csv")
write_csv(res$peaks, "data_inter/disproportionate_cohorts.csv")
write_csv(res$pairs, "data_inter/block_pairs.csv")
write_csv(tab_s5_1, "data_inter/tableS5_1_exceptional_by_age_sex.csv")

# alternative reference set (Supplementary Material S1) ------------------------
load("data_inter/cohort_conts_gap_changes_altref.RData")

res_alt <- exceptionality_criterion(dt_gaps2_alt, k = k_each_side, z = z_thr,
                                    min_cts = min_cts, min_ags = min_ags, L_max = L_max)

cat("\n--- alternative reference set (", paste(unique(dt_gaps2_alt$cts), collapse = ", "), ")\n",
    "tested contributions:", res_alt$n_tested, " exceptional:", res_alt$n_exceptional,
    sprintf("(%.1f%%)", res_alt$pct_exceptional), "\n",
    "replicated single peaks:", format_peaks(res_alt$peaks), "\n",
    "opposite-signed pairs found:", nrow(res_alt$pairs),
    "; replicated pairs:", format_blocks(res_alt$blocks), "\n")

write_csv(res_alt$exceptional, "data_inter/exceptional_cohorts_altref.csv")

# power check on the simulated scenarios ---------------------------------------
# One simulated comparator, so only the replication across age intervals
# applies. Scenario B (cohorts 1945-1954 with 15% lower mortality) should give
# a narrowing peak at 1944 (the cohort replaced by the first advantaged one)
# and a widening peak at 1954 (the last advantaged cohort) in every age
# interval where the block is observed, and the pair (1944, 1954).
sim <-
  read_rds("data_inter/simulation_scc.rds") %>%
  mutate(scen = str_sub(as.character(scen), 1, 1), cts = "SIM", sex = "t")

flags_sim <- flag_exceptional(sim, c("scen", "ags"), k = k_each_side, z = z_thr) %>%
  filter(exceptional) %>%
  select(scen, ags, cohort, scc, z_score, direction)

cat("\n--- simulation: exceptional cohorts by scenario and age interval\n")
print(flags_sim %>%
        group_by(scen, ags) %>%
        summarise(flagged = paste0(cohort, ifelse(direction == "narrowing", "(-)", "(+)"), collapse = " "),
                  .groups = "drop"), n = 40)

cat("\n--- simulation: single peaks recurring in >=", min_ags, "age intervals\n")
print(flags_sim %>% reframe(n_ags = n_distinct(ags), .by = c(scen, cohort, direction)) %>%
        filter(n_ags >= min_ags) %>% arrange(scen, cohort))

pairs_sim <- block_pairs(flags_sim %>% mutate(exceptional = TRUE, cts = "SIM", sex = "t"),
                         c("scen", "ags"), L_max)
cat("\n--- simulation: block pairs recurring in >=", min_ags, "age intervals\n")
print(pairs_sim %>% reframe(n_ags = n_distinct(ags), .by = c(scen, block, entry, exit)) %>%
        filter(n_ags >= min_ags))

write_csv(flags_sim, "data_inter/exceptional_cohorts_simulation.csv")

# sensitivity (Table S5.3) -------------------------------------------------------
settings <- tribble(
  ~scale, ~k, ~z,
  "mad",  5,  3,      # main
  "mad",  5,  2.5,
  "mad",  5,  3.5,
  "mad",  8,  3,
  "sd",   5,  2,      # first version of the criterion
  "sd",   5,  2.5,
  "sd",   8,  2
)

sens <- map_dfr(seq_len(nrow(settings)), function(i){
  s <- settings[i, ]
  r <- exceptionality_criterion(dt_gaps2, k = s$k, z = s$z, scale = s$scale,
                                min_cts = min_cts, min_ags = min_ags, L_max = L_max)
  tibble(scale = s$scale, k_each_side = s$k, z_thr = s$z,
         pct_exceptional = round(r$pct_exceptional, 1),
         peaks = format_peaks(r$peaks),
         blocks = format_blocks(r$blocks))
})

cat("\n--- sensitivity to the scale estimator, the window and the threshold\n")
print(as.data.frame(sens))

write_csv(sens, "data_inter/exceptional_cohorts_sensitivity.csv")
