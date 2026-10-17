# ==============================================================================
# 11. Formal criterion for exceptional cohort contributions (Methods, Step 3).
#
# A cohort contribution (unscaled, dt_coh) is "exceptional" when it deviates
# from the mean of its neighbouring cohorts in the same age group and
# comparator by more than z standard deviations of those neighbours. A cohort
# pattern is "disproportionate" when, for the same sex and cohort, this occurs
# (with the same sign) in at least min_cts comparator countries in at least
# min_ags age groups. The parameters below are the ones in the manuscript
# (shown there between brackets).
#
# Input : data_inter/cohort_conts_gap_changes.RData (07_fig5_cohort_contributions_gap.R)
# Output: data_inter/exceptional_cohorts.csv (all flagged contributions)
#         data_inter/disproportionate_cohorts.csv (cohorts meeting the criterion)
# ==============================================================================
rm(list = ls()); gc()
library(tidyverse)

load("data_inter/cohort_conts_gap_changes.RData")

k_each_side <- 5   # 10 neighbouring cohorts: 5 on each side
z           <- 2   # standard deviations
min_ags     <- 2   # age groups
min_cts     <- 3   # comparator countries (out of 4)

flags <-
  dt_gaps2 %>%
  filter(sex != "t", ags != "100_109") %>%
  group_by(cts, sex, ags) %>%
  arrange(cohort, .by_group = TRUE) %>%
  mutate(nb_mean = sapply(row_number(), function(i){
           j <- setdiff(max(1, i - k_each_side):min(n(), i + k_each_side), i)
           mean(dt_coh[j]) }),
         nb_sd = sapply(row_number(), function(i){
           j <- setdiff(max(1, i - k_each_side):min(n(), i + k_each_side), i)
           sd(dt_coh[j]) }),
         dev = dt_coh - nb_mean,
         exceptional = abs(dev) > z * nb_sd,
         direction = ifelse(dev < 0, "narrowing", "widening")) %>%
  ungroup() %>%
  select(cts, sex, ags, cohort, dt_coh, scc, nb_mean, nb_sd, dev, exceptional, direction)

exceptional <- flags %>% filter(exceptional)

cat("exceptional contributions:", nrow(exceptional), "of", nrow(flags),
    sprintf("(%.1f%%)", 100 * nrow(exceptional) / nrow(flags)), "\n")

# same sex, cohort, age group and direction in >= min_cts comparators ...
by_ags <-
  exceptional %>%
  reframe(n_cts = n_distinct(cts), .by = c(sex, cohort, ags, direction)) %>%
  filter(n_cts >= min_cts)

# ... and that in >= min_ags age groups
disprop <-
  by_ags %>%
  reframe(n_ags = n_distinct(ags),
          age_groups = paste(sort(as.character(ags)), collapse = ", "),
          .by = c(sex, cohort, direction)) %>%
  filter(n_ags >= min_ags) %>%
  arrange(sex, cohort)

cat("\ncohort x sex x direction meeting the criterion (>=", min_cts, "comparators in >=", min_ags, "age groups):\n")
print(as.data.frame(disprop))

cat("\nexceptional contributions by age group and sex (any comparator):\n")
print(exceptional %>% count(sex, ags) %>% pivot_wider(names_from = sex, values_from = n, values_fill = 0))

write_csv(exceptional, "data_inter/exceptional_cohorts.csv")
write_csv(disprop, "data_inter/disproportionate_cohorts.csv")

# ------------------------------------------------------------------------------
# Power check on the simulated scenarios of Figure S3 (08_figS3_simulation_scc.R)
#
# The same two-sided criterion is applied to the simulated cohort contributions.
# In scenario B (cohorts 1945-1954 with 15% lower mortality) it should flag the
# two edges of the block: the cohort replaced by the first advantaged one (1944,
# narrowing) and the last advantaged cohort (1954, widening). In the period-only
# scenarios (A and C) any flagged cohort is an edge of the period band, which
# falls on a different cohort in each age interval and therefore never meets
# the replication requirement across age intervals (one simulated comparator
# only, so the comparator condition cannot be applied).
# Output: data_inter/exceptional_cohorts_simulation.csv
# ------------------------------------------------------------------------------
sim <- 
  read_rds("data_inter/simulation_scc.rds") %>% 
  mutate(scen = str_sub(as.character(scen), 1, 1))

flags_sim <- 
  sim %>% 
  group_by(scen, ags) %>% 
  arrange(cohort, .by_group = TRUE) %>% 
  mutate(nb_mean = sapply(row_number(), function(i){
           j <- setdiff(max(1, i - k_each_side):min(n(), i + k_each_side), i)
           mean(dt_coh[j]) }),
         nb_sd = sapply(row_number(), function(i){
           j <- setdiff(max(1, i - k_each_side):min(n(), i + k_each_side), i)
           sd(dt_coh[j]) }),
         dev = dt_coh - nb_mean,
         exceptional = abs(dev) > z * nb_sd,
         direction = ifelse(dev < 0, "narrowing", "widening")) %>% 
  ungroup() %>% 
  filter(exceptional) %>% 
  select(scen, ags, cohort, scc, dev, nb_sd, direction)

cat("\n--- simulation: exceptional cohorts by scenario and age interval\n")
print(flags_sim %>% 
        group_by(scen, ags) %>% 
        summarise(flagged = paste0(cohort, ifelse(direction == "narrowing", "(-)", "(+)"), collapse = " "), 
                  .groups = "drop"), n = 40)

cat("\n--- simulation: cohorts exceptional in >=", min_ags, "age intervals (same direction)\n")
print(flags_sim %>% 
        reframe(n_ags = n_distinct(ags), .by = c(scen, cohort, direction)) %>% 
        filter(n_ags >= min_ags) %>% 
        arrange(scen, cohort))

write_csv(flags_sim, "data_inter/exceptional_cohorts_simulation.csv")
