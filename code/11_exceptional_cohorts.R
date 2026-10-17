# ==============================================================================
# 11. Operational criterion for exceptional and disproportionate cohort
#     contributions (Methods, Step 3), its power check on the simulation of
#     Figure S3, and a sensitivity analysis.
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
# Definitions (parameters below).
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
#         data_inter/simulation_scc.rds (08_figS3_simulation_scc.R)
# Output: data_inter/exceptional_cohorts.csv            all flagged contributions
#         data_inter/disproportionate_cohorts.csv       single peaks meeting the criterion
#         data_inter/block_pairs.csv                    opposite-signed pairs found
#         data_inter/exceptional_cohorts_simulation.csv power check on the simulation
#         data_inter/exceptional_cohorts_sensitivity.csv results for other settings
# ==============================================================================
rm(list = ls()); gc()
library(tidyverse)

# parameters -------------------------------------------------------------------
k_each_side <- 5    # neighbouring cohorts on each side (10 in total)
z_thr       <- 3    # threshold on the robust z-score
min_ags     <- 2    # age intervals in which the signal must recur
min_cts     <- 3    # comparator countries in which it must recur (out of 4)
L_max       <- 15   # longest block of cohorts considered for the pair signature

# functions --------------------------------------------------------------------
# robust z-score of each value against its neighbours (excluding itself);
# NA for the first and last k cohorts (one-sided neighbourhoods are not tested)
robust_z <- function(x, k = k_each_side){
  n <- length(x)
  sapply(seq_len(n), function(i){
    if (i <= k | i > n - k) return(NA_real_)
    j <- setdiff((i - k):(i + k), i)
    (x[i] - median(x[j])) / (1.4826 * mad(x[j], constant = 1))
  })
}

# classical version used in a first version of the criterion (sensitivity)
sd_z <- function(x, k = k_each_side){
  n <- length(x)
  sapply(seq_len(n), function(i){
    if (i <= k | i > n - k) return(NA_real_)
    j <- setdiff((i - k):(i + k), i)
    (x[i] - mean(x[j])) / sd(x[j])
  })
}

# flag exceptional contributions within each group
flag_exceptional <- function(dt, groups, k = k_each_side, z = z_thr, scale = c("mad", "sd")){
  scale <- match.arg(scale)
  dt %>%
    group_by(across(all_of(groups))) %>%
    arrange(cohort, .by_group = TRUE) %>%
    mutate(z_score = if (scale == "mad") robust_z(dt_coh, k) else sd_z(dt_coh, k),
           exceptional = !is.na(z_score) & abs(z_score) > z,
           direction = ifelse(z_score < 0, "narrowing", "widening")) %>%
    ungroup()
}

# single peaks recurring across comparators and age intervals
disproportionate_peaks <- function(flags, min_cts, min_ags){
  flags %>%
    filter(exceptional) %>%
    reframe(n_cts = n_distinct(cts), .by = c(sex, cohort, ags, direction)) %>%
    filter(n_cts >= min_cts) %>%
    reframe(n_ags = n_distinct(ags),
            age_groups = paste(sort(as.character(ags)), collapse = ", "),
            .by = c(sex, cohort, direction)) %>%
    filter(n_ags >= min_ags) %>%
    arrange(sex, cohort)
}

# opposite-signed pairs within each group: entry and exit of a block
block_pairs <- function(flags, groups, L_max){
  flags %>%
    filter(exceptional) %>%
    group_by(across(all_of(groups))) %>%
    group_modify(function(d, g){
      a <- d$cohort[d$direction == "narrowing"]
      b <- d$cohort[d$direction == "widening"]
      adv <- expand_grid(entry = a, exit = b) %>% filter(exit - entry >= 2, exit - entry <= L_max) %>% mutate(block = "advantaged")
      dis <- expand_grid(entry = b, exit = a) %>% filter(exit - entry >= 2, exit - entry <= L_max) %>% mutate(block = "disadvantaged")
      bind_rows(adv, dis)
    }) %>%
    ungroup() %>%
    mutate(length = exit - entry)
}

# pairs recurring across comparators and age intervals
disproportionate_blocks <- function(pairs, min_cts, min_ags){
  pairs %>%
    reframe(n_cts = n_distinct(cts), .by = c(sex, block, entry, exit, ags)) %>%
    filter(n_cts >= min_cts) %>%
    reframe(n_ags = n_distinct(ags),
            age_groups = paste(sort(as.character(ags)), collapse = ", "),
            .by = c(sex, block, entry, exit)) %>%
    filter(n_ags >= min_ags) %>%
    arrange(sex, entry)
}

# observed data ----------------------------------------------------------------
load("data_inter/cohort_conts_gap_changes.RData")

base <- dt_gaps2 %>% filter(sex != "t", ags != "100_109")
grp  <- c("cts", "sex", "ags")

flags <- flag_exceptional(base, grp) %>%
  select(cts, sex, ags, cohort, dt_coh, scc, z_score, exceptional, direction)

tested <- flags %>% filter(!is.na(z_score))
exceptional <- flags %>% filter(exceptional)

cat("tested contributions:", nrow(tested), " exceptional:", nrow(exceptional),
    sprintf("(%.1f%%)", 100 * nrow(exceptional) / nrow(tested)), "\n")

disprop <- disproportionate_peaks(flags, min_cts, min_ags)
cat("\nsingle peaks meeting the criterion (same sign in >=", min_cts, "comparators and >=", min_ags, "age intervals):\n")
print(as.data.frame(disprop))

pairs <- block_pairs(flags, grp, L_max)
blocks <- disproportionate_blocks(pairs, min_cts, min_ags)
cat("\nopposite-signed pairs found in the data:", nrow(pairs),
    "; pairs meeting the replication requirement:", nrow(blocks), "\n")
if (nrow(blocks) > 0) print(as.data.frame(blocks))

cat("\nexceptional contributions by age interval and sex (any comparator):\n")
print(exceptional %>% count(sex, ags) %>% pivot_wider(names_from = sex, values_from = n, values_fill = 0))

write_csv(exceptional, "data_inter/exceptional_cohorts.csv")
write_csv(disprop, "data_inter/disproportionate_cohorts.csv")
write_csv(pairs, "data_inter/block_pairs.csv")

# power check on the simulated scenarios ---------------------------------------
# One simulated comparator, so only the replication across age intervals
# applies. Scenario B (cohorts 1945-1954 with 15% lower mortality) should give
# a narrowing peak at 1944 (the cohort replaced by the first advantaged one)
# and a widening peak at 1954 (the last advantaged cohort) in every age
# interval where the block is observed, and the pair (1944, 1954).
sim <-
  read_rds("data_inter/simulation_scc.rds") %>%
  mutate(scen = str_sub(as.character(scen), 1, 1), cts = "SIM", sex = "t")

flags_sim <- flag_exceptional(sim, c("scen", "ags")) %>%
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

# sensitivity --------------------------------------------------------------------
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
  f <- flag_exceptional(base, grp, k = s$k, z = s$z, scale = s$scale)
  d <- disproportionate_peaks(f, min_cts, min_ags)
  b <- disproportionate_blocks(block_pairs(f, grp, L_max), min_cts, min_ags)
  tibble(scale = s$scale, k_each_side = s$k, z_thr = s$z,
         pct_exceptional = round(100 * mean(f$exceptional[!is.na(f$z_score)]), 1),
         peaks = ifelse(nrow(d) == 0, "none", paste(sprintf("%s%d(%s)", d$sex, d$cohort, d$age_groups), collapse = "; ")),
         blocks = ifelse(nrow(b) == 0, "none", paste(sprintf("%s %d-%d", b$sex, b$entry, b$exit), collapse = "; ")))
})

cat("\n--- sensitivity to the scale estimator, the window and the threshold\n")
print(as.data.frame(sens))

write_csv(sens, "data_inter/exceptional_cohorts_sensitivity.csv")
