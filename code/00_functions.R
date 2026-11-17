# ==============================================================================
# Functions shared by the manuscript scripts
# "From post-war hardship to longevity leader: a period and cohort analysis of
#  the Spanish mortality paradox"
# ==============================================================================
library(tidyverse)
options(scipen = 999)

# ------------------------------------------------------------------------------
# Typeface of the figures: a serif family (Times-like, as in the manuscript),
# set once here and picked up by every script that sources this file. Each
# script adds a complete theme (theme_bw() or theme_minimal()) to its plots,
# which would reset the family to the device default, so both themes are
# redefined with a serif default; geom_text() and annotate("text") inherit the
# family from the theme (ggplot2 >= 4.0).
fig_family <- "serif"
theme_set(theme_bw(base_family = fig_family))
theme_bw <- function(..., base_family = fig_family) ggplot2::theme_bw(..., base_family = base_family)
theme_minimal <- function(..., base_family = fig_family) ggplot2::theme_minimal(..., base_family = base_family)

# ------------------------------------------------------------------------------
# Arriaga (1984) decomposition of the difference between two (partial) life
# expectancies into single-age contributions, population 1 taken as the
# starting point. Internal helper used by est_age_cont().
#
# lt1, lt2: life tables restricted to the age window, with lx, Lx, Tx and ex
#           recomputed within the window (see est_age_cont()).
#   dd : direct effect (change in years lived within the age interval)
#   di : indirect + interaction effect (change in years lived beyond the
#        interval because of survivors added to the next age)
#   dt : total contribution of the age to the difference ex2 - ex1
# ------------------------------------------------------------------------------
arriaga_dir <- function(lt1, lt2, amin, amax){

  # radix at the first age of the window
  l01 <-
    lt1 %>%
    filter(age == amin) %>%
    pull(lx)

  lt1 %>%
    select(age, lx1 = lx, Lx1 = Lx, Tx1 = Tx, ex1 = ex) %>%
    left_join(lt2 %>%
                select(age, lx2 = lx, Lx2 = Lx, Tx2 = Tx, ex2 = ex),
              by = join_by(age)) %>%
    mutate(dd = ifelse(lx1 >0 & lx2 >0,
                       lx1/l01*(Lx2/lx2-Lx1/lx1),
                       0),
           di = ifelse(age < amax & lx2 > 0 & lead(lx2) > 0,
                       lead(Tx2)/l01*(lx1/lx2-lead(lx1)/lead(lx2)),
                       0),
           dt = ifelse(age < amax, dd + di, dd))
}

# ------------------------------------------------------------------------------
# Age-specific contributions to the difference e2 - e1 between the (partial)
# life expectancies of population 1 (c1, s1, y1) and population 2 (c2, s2, y2).
#
# Life tables are taken from `hmd_lt` (HMD complete life tables, see
# 01_hmd_life_tables.R), restricted to ages amin:amax. Tx and ex are recomputed
# within the age window, so e_amin is the partial life expectancy between
# ages amin and amax + 1. The contributions (dt) sum to e2 - e1.
#
# symmetric = TRUE (default, used in the manuscript): Arriaga's decomposition
# is asymmetric (it depends on which population is the starting point), so the
# contributions are averaged over the two directions (Andreev et al. 2002;
# Ponnapalli 2005): (contributions to e2 - e1 starting from population 1
# minus contributions to e1 - e2 starting from population 2) / 2.
# symmetric = FALSE returns the plain Arriaga decomposition from population 1.
# ------------------------------------------------------------------------------
est_age_cont <- function(c1 = "ESP", s1 = "t", y1 = 1950,
                         c2 = "ESP", s2 = "t", y2 = 2000,
                         amin = 10, amax = 85,
                         symmetric = TRUE){

  lt1 <-
    hmd_lt %>%
    filter(code == c1,
           year == y1,
           age %in% amin:amax,
           sex == s1) %>%
    mutate(Tx = order_by(-age, cumsum(Lx)),
           ex = ifelse(lx>0, round(Tx/lx, 2), 0))

  lt2 <-
    hmd_lt %>%
    filter(code == c2,
           year == y2,
           age %in% amin:amax,
           sex == s2) %>%
    mutate(Tx = order_by(-age, cumsum(Lx)),
           ex = ifelse(lx>0, round(Tx/lx, 2), 0))

  # contributions to e2 - e1, population 1 as starting point
  fwd <- arriaga_dir(lt1, lt2, amin, amax)

  if (!symmetric) return(fwd)

  # contributions to e1 - e2, population 2 as starting point
  bwd <- arriaga_dir(lt2, lt1, amin, amax)

  fwd %>%
    mutate(dd = (dd - bwd$dd) / 2,
           di = (di - bwd$di) / 2,
           dt = (dt - bwd$dt) / 2)
}

# ------------------------------------------------------------------------------
# Cohort-specific contributions to the change in the (partial) life expectancy
# gap between two countries (c2 - c1, i.e. REF - ESP when c1 = "ESP") over the
# period ymin:ymax. Negative values = the gap narrows.
#
# For each year the gap is decomposed by age with est_age_cont(). The change in
# the contribution of age x between consecutive years y and y + 1 is assigned
# to cohort c = y - x (the cohort aged x in year y), and the changes are then
# cumulated along the Lexis diagonal of each cohort (manuscript, Step 3):
#     Cc_{c,[a,a+n)} = sum_{t=0}^{n-1} [ D_{a+t}(c+a+t+1) - D_{a+t}(c+a+t) ]
# Only cohorts observed at every age of the window, with a defined change at
# each of them, are kept (parallelograms in Figure 2 of the manuscript), i.e.
# cohorts ymin - amin <= c <= ymax - amax - 1.
#
# Returns a list:
#   long     : age x year x sex contributions and their annual changes
#   coh_cont : one row per cohort and sex (dt_coh = cumulated contribution to
#              the change in the gap; negative = gap narrows)
#   summ     : partial life expectancies of both countries by year and sex
# ------------------------------------------------------------------------------
est_coh_cont_diff_countries <- function(c1 = "ESP", c2 = "SWE",
                                        amin = 10, amax = 85,
                                        ymin = 1930, ymax = 2023,
                                        sexes = c("m", "f", "t"),
                                        symmetric = TRUE){

  diffs <- tibble()
  summs <- tibble()
  for(i in ymin:ymax){

    for (s in sexes){

      tmp1 <-
        est_age_cont(c1 = c1, s1 = s, y1 = i,
                     c2 = c2, s2 = s, y2 = i,
                     amin = amin, amax = amax,
                     symmetric = symmetric)

      tmp2 <-
        tmp1 %>%
        select(age, dd, di, dt) %>%
        mutate(sex = s,
               year = i)

      diffs <-
        bind_rows(diffs, tmp2)

      summ1 <-
        tmp1 %>%
        mutate(ags = paste0(min(age), "_", max(age))) %>%
        filter(age == min(age)) %>%
        select(ags, ex1, ex2) %>%
        mutate(year = i,
               sex = s,
               cts = paste0(c1, "_", c2))

      summs <-
        summs %>%
        bind_rows(summ1)

    }
  }

  dt <-
    diffs %>%
    group_by(sex, age) %>%
    arrange(sex, age, year) %>%
    # change in the age-specific contribution from this year to the next,
    # assigned to the cohort aged `age` in `year`
    mutate(dt_chg = lead(dt) - dt) %>%
    ungroup() %>%
    filter(!is.na(dt_chg)) %>%
    mutate(cohort = year - age) %>%
    group_by(sex) %>%
    mutate(dt_tot_sqr = sum(dt_chg)) %>%
    # keep only cohorts with a change defined at all ages of the window
    group_by(cohort, sex) %>%
    mutate(yrs = n()) %>%
    ungroup() %>%
    filter(yrs == amax - amin + 1) %>%
    arrange(sex, cohort) %>%
    group_by(sex) %>%
    mutate(pos = ifelse(dt_chg > 0 , dt_chg, 0), # contributions widening the gap
           neg = ifelse(dt_chg < 0 , dt_chg, 0), # contributions narrowing the gap
           dt_tot = sum(dt_chg),
           pos_tot = sum(pos),
           neg_tot = sum(neg),
           neg_rel = 100*neg/neg_tot) %>% # relative contribution to the narrowing
    ungroup() %>%
    mutate(yrs = paste0(ymin, "_", ymax),
           ags = paste0(amin, "_", amax),
           chs = paste0(min(cohort), "_", max(cohort)),
           cts = paste0(c1, "_", c2))

  coh_cont <-
    dt %>%
    reframe(dt_coh = sum(dt_chg),
            dt_pos = sum(pos),
            dt_neg = sum(neg),
            dt_tot = mean(dt_tot),
            dt_pos_tot = mean(pos_tot),
            dt_neg_tot = mean(neg_tot),
            dt_rel = sum(neg_rel),
            dt_tot_sqr = mean(dt_tot_sqr),
            .by = c(cohort, sex)) %>%
    arrange(sex, cohort) %>%
    mutate(yrs = paste0(ymin, "_", ymax),
           ags = paste0(amin, "_", amax),
           chs = paste0(min(cohort), "_", max(cohort)),
           cts = paste0(c1, "_", c2))

  out <- list(long = dt, coh_cont = coh_cont, summ = summs)
  return(out)

}

# ------------------------------------------------------------------------------
# Decomposition of the change in the partial life expectancy gap between
# country c1 and each reference country in cts_ref, by age groups, over
# ymin:ymax (Figures 5 and S1). Wraps est_coh_cont_diff_countries().
# ymax can be a single year or a named vector with one year per reference
# country (e.g. c(DNK = 2023, GBR_NP = 2022) when a series ends earlier).
#
# Returns a list:
#   cohs : one row per cohort x sex x age group x comparator, with the
#          cumulated contribution (dt_coh), its scaled version (scc, divided
#          by the largest absolute contribution within sex x age group x
#          comparator) and per2 (= cohort when the cohort crossed the age
#          group at some point during mark_years, NA otherwise)
#   summ : partial life expectancies of both countries by year and sex
# ------------------------------------------------------------------------------
decompose_gap_by_age_groups <- function(c1 = "ESP",
                                        cts_ref = c("SWE", "CHE", "NOR", "NLD"),
                                        ymin = 1940, ymax = 2023,
                                        age_starts = seq(0, 100, 10), width = 10,
                                        sexes = c("m", "f", "t"),
                                        mark_years = 1940:1960){

  cohs <- tibble()
  summ <- tibble()
  for(a in age_starts){
    for(ct in cts_ref){
      # last year of this comparison (ymax may be a named vector, see above)
      ymax_ct <- if (is.null(names(ymax))) ymax else ymax[[ct]]
      cat(paste0("ages ", a, "-", a + width - 1, " ", c1, " vs ", ct,
                 " (", ymin, "-", ymax_ct, ")\n"))
      res <-
        est_coh_cont_diff_countries(c1 = c1, c2 = ct,
                                    amin = a, amax = a + width - 1,
                                    ymin = ymin, ymax = ymax_ct,
                                    sexes = sexes)
      summ <- bind_rows(summ, res$summ)
      cohs <- bind_rows(cohs, res$coh_cont)
    }
  }

  ag_levels <- paste0(age_starts, "_", age_starts + width - 1)

  cohs <-
    cohs %>%
    mutate(ags = factor(ags, levels = ag_levels),
           amin = as.integer(str_extract(as.character(ags), "^[0-9]+")),
           # cohorts that crossed the age group at some point during mark_years
           per2 = ifelse(cohort + amin <= max(mark_years) &
                           cohort + amin + width - 1 >= min(mark_years),
                         cohort, NA)) %>%
    group_by(sex, ags, cts) %>%
    mutate(scc = dt_coh / max(abs(dt_coh))) %>%
    ungroup()

  list(cohs = cohs, summ = summ)
}

# ------------------------------------------------------------------------------
# Figure 5 / S1: scaled cohort contributions by sex and age group, one colour
# per comparator country; dark marks at the top of each panel for the cohorts
# that crossed the age group during the marked period (per2). `labels` (a named
# vector, HMD code -> name) replaces the codes in the legend.
# ------------------------------------------------------------------------------
plot_scc <- function(cohs,
                     cols = c("#AE2012", "#94D2BD", "#005F73", "#EE9B00"),
                     legend_title = "Gap between Spain and:",
                     drop_ags = "100_109",
                     labels = NULL){

  cohs %>%
    filter(sex != "t", !(ags %in% drop_ags)) %>%
    mutate(sex = ifelse(sex == "f", "Female", "Male"),
           ags = str_replace(as.character(ags), "_", "-"),
           ags = fct_reorder(ags, amin),
           cts = str_remove(cts, "^[A-Z_]+?_"),
           cts = if (is.null(labels)) cts else factor(unname(labels[cts]), levels = labels)) %>%
    ggplot(aes(cohort, scc))+
    facet_grid(fct_rev(ags)~sex, scales = "free", space = "free", switch = "y")+
    geom_line(aes(cohort, scc, col = cts), alpha = .3, lwd = 0.1)+
    geom_point(aes(cohort, scc, col = cts), alpha = .5, size = 0.03)+
    geom_point(aes(x = per2, y = Inf), col = "#033b86", size = .5, na.rm = TRUE)+
    geom_smooth(aes(group = 1), method = "loess", span = 0.4, se = FALSE,
                lwd = 0.3, alpha = .6, color = "black") +
    geom_hline(yintercept = 0, lty = "dashed", lwd = .2)+
    scale_x_continuous(breaks = seq(1800, 2020, 10))+
    scale_y_continuous(position = "right", breaks = c(-1, 0, 1))+
    scale_color_manual(values = cols)+
    labs(y = "Scaled cohort contribution (SCC) to the LE gap change",
         x = "Cohort",
         color = legend_title)+
    coord_cartesian(ylim = c(-1,1), expand = 0)+
    guides(color = guide_legend(override.aes = list(size = 2.5, linewidth = 1.5)))+
    theme_bw()+
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 9),
          panel.spacing.x = unit(0.5, "cm"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
          axis.title.y.right = element_text(vjust = 0.5, angle = 90),
          legend.position = "bottom",
          legend.key.size = unit(.8, "cm"))
}

# ==============================================================================
# Detrended cohort effects (Figures 3, S2 and S8.1; Tables S2.1 and S8.1)
# ==============================================================================

# ------------------------------------------------------------------------------
# Age-drift-Period-Cohort model (Epi::apc.fit, factor model, "AdPC"
# parametrization, drift extracted with the unweighted inner product) fitted to
# single-year death counts and exposures (data_inter/deaths_exposures_hmd.rds,
# see 02_hmd_deaths_exposures.R) for one country, sex, age range and period.
# Returns the age, period and cohort effects and the drift.
# ------------------------------------------------------------------------------
est_apc <- function(hmd, cd = "ESP", sx = "female",
                    amin = 0, amax = 99,
                    pmin = 1940, pmax = 2023) {
  dt <-
    hmd %>%
    filter(code == cd, sex == sx, age %in% amin:amax, year %in% pmin:pmax) %>%
    rename(D = dts, Y = exposure, A = age, P = year) %>%
    mutate(C = P - A)

  # Poisson APC model with age, period and cohort as factors; drift extracted
  # with the unweighted (unit) inner product; period and cohort effects
  # detrended (zero mean and zero slope)
  apc_c <-
    Epi::apc.fit(dt, model = "factor", dr.extr = "1", parm = "AdPC", scale = 10^5,
                 print.AOV = FALSE)

  list(age = apc_c$Age %>% as_tibble() %>% mutate(code = cd, sex = sx),
       per = apc_c$Per %>% as_tibble() %>% mutate(code = cd, sex = sx),
       coh = apc_c$Coh %>% as_tibble() %>% mutate(code = cd, sex = sx),
       drift = apc_c$Drift[1, 1])
}

# ------------------------------------------------------------------------------
# Detrended cohort effects for several countries, age groups and both sexes,
# one row per cohort (C-RR = rate ratio, lrr = its log), as saved in
# data_inter/apc_cohort_effects.rds. One factor model per country x age group
# x sex (30 models for five countries; 15-20 minutes).
# ------------------------------------------------------------------------------
fit_cohort_effects <- function(hmd, cts, names_cts,
                               age_groups = list(c(0, 49), c(50, 99), c(0, 99)),
                               pmin = 1940, pmax = 2023) {
  cohs <- tibble()
  for (cd in cts) {
    for (ag in age_groups) {
      for (sx in c("female", "male")) {
        cat(paste0(cd, " ", sx, " ages ", ag[1], "-", ag[2], " (", pmin, "-", pmax, ")\n"))
        res <- est_apc(hmd, cd = cd, sx = sx, amin = ag[1], amax = ag[2],
                       pmin = pmin, pmax = pmax)
        cohs <-
          bind_rows(cohs,
                    res$coh %>% mutate(ages = paste0(ag[1], "-", ag[2]),
                                       drift = res$drift))
      }
    }
  }

  cohs %>%
    mutate(ages = factor(ages, levels = sapply(age_groups, function(ag) paste0(ag[1], "-", ag[2]))),
           sex = str_to_title(sex),
           country = factor(names_cts[code], levels = names_cts),
           lrr = log(`C-RR`),
           period = paste0(pmin, "-", pmax))
}

# ------------------------------------------------------------------------------
# Local extrema of the detrended cohort curves (Table S2.1). Because the drift
# extraction only tilts the detrended curve linearly, the depth of an extremum
# is measured relative to the straight line joining the neighbouring turning
# points of the curve (prominence, in log RR), which is invariant to that
# tilt. The deviation from the mean of the +/-10 neighbouring cohorts is also
# reported.
# ------------------------------------------------------------------------------
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

# features reported in Table S2.1: the 1940s minimum and the 1960s maximum
extrema_features <- tribble(
  ~feature,            ~type, ~window,
  "1940s minimum",     "min", 1935:1955,
  "1960s maximum",     "max", 1955:1970
)

# Table S2.1: location and depth of each feature by country, sex and age group
cohort_extrema_table <- function(cohs, features = extrema_features){
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
}

# ==============================================================================
# Criterion for exceptional and disproportionate cohort contributions
# (Methods, Step 3; script 11; Tables S5.1-S5.3 and S8.1)
# ==============================================================================

# ------------------------------------------------------------------------------
# Robust z-score of each value against its k neighbours on each side (excluding
# itself): (x - median of neighbours) / (1.4826 * MAD of neighbours). NA for
# the first and last k cohorts (one-sided neighbourhoods are not tested).
# ------------------------------------------------------------------------------
robust_z <- function(x, k = 5){
  n <- length(x)
  sapply(seq_len(n), function(i){
    if (i <= k | i > n - k) return(NA_real_)
    j <- setdiff((i - k):(i + k), i)
    (x[i] - median(x[j])) / (1.4826 * mad(x[j], constant = 1))
  })
}

# classical version (mean and standard deviation of the neighbours) used in a
# first version of the criterion; kept for the sensitivity analysis
sd_z <- function(x, k = 5){
  n <- length(x)
  sapply(seq_len(n), function(i){
    if (i <= k | i > n - k) return(NA_real_)
    j <- setdiff((i - k):(i + k), i)
    (x[i] - mean(x[j])) / sd(x[j])
  })
}

# flag exceptional contributions (|z| > z) within each group
flag_exceptional <- function(dt, groups, k = 5, z = 3, scale = c("mad", "sd")){
  scale <- match.arg(scale)
  dt %>%
    group_by(across(all_of(groups))) %>%
    arrange(cohort, .by_group = TRUE) %>%
    mutate(z_score = if (scale == "mad") robust_z(dt_coh, k) else sd_z(dt_coh, k),
           exceptional = !is.na(z_score) & abs(z_score) > z,
           direction = ifelse(z_score < 0, "narrowing", "widening")) %>%
    ungroup()
}

# single peaks recurring (same cohort, sex and sign) in at least min_cts
# comparators within an age interval and in at least min_ags age intervals
disproportionate_peaks <- function(flags, min_cts = 3, min_ags = 2){
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

# opposite-signed pairs within each group: a narrowing peak at cohort c and a
# widening peak at c + L, 2 <= L <= L_max (entry and exit of an advantaged
# block), or the reverse (a disadvantaged block)
block_pairs <- function(flags, groups, L_max = 15){
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
disproportionate_blocks <- function(pairs, min_cts = 3, min_ags = 2){
  pairs %>%
    reframe(n_cts = n_distinct(cts), .by = c(sex, block, entry, exit, ags)) %>%
    filter(n_cts >= min_cts) %>%
    reframe(n_ags = n_distinct(ags),
            age_groups = paste(sort(as.character(ags)), collapse = ", "),
            .by = c(sex, block, entry, exit)) %>%
    filter(n_ags >= min_ags) %>%
    arrange(sex, entry)
}

# ------------------------------------------------------------------------------
# The whole criterion applied to the cohort contributions of
# decompose_gap_by_age_groups() (one sex x age interval x comparator series
# each; both sexes combined and the open age group are not tested).
# Returns the flagged contributions, the single peaks and pairs meeting the
# replication requirement, and the counts reported in Tables S5.3 and S8.1.
# ------------------------------------------------------------------------------
exceptionality_criterion <- function(cohs, k = 5, z = 3, scale = "mad",
                                     min_cts = 3, min_ags = 2, L_max = 15,
                                     drop_ags = "100_109"){
  base <- cohs %>% filter(sex != "t", !(ags %in% drop_ags))
  grp  <- c("cts", "sex", "ags")
  flags <-
    flag_exceptional(base, grp, k = k, z = z, scale = scale) %>%
    select(cts, sex, ags, cohort, dt_coh, scc, z_score, exceptional, direction)
  pairs <- block_pairs(flags, grp, L_max)
  tested <- flags$exceptional[!is.na(flags$z_score)]
  list(flags = flags,
       exceptional = flags %>% filter(exceptional),
       peaks = disproportionate_peaks(flags, min_cts, min_ags),
       pairs = pairs,
       blocks = disproportionate_blocks(pairs, min_cts, min_ags),
       n_tested = length(tested),
       n_exceptional = sum(tested),
       pct_exceptional = 100 * mean(tested))
}

# compact text versions of the replicated peaks and pairs (tables)
format_peaks <- function(d){
  if (nrow(d) == 0) return("none")
  paste(sprintf("%s%d(%s)", d$sex, d$cohort, d$age_groups), collapse = "; ")
}
format_blocks <- function(b){
  if (nrow(b) == 0) return("none")
  paste(sprintf("%s %d-%d", b$sex, b$entry, b$exit), collapse = "; ")
}
