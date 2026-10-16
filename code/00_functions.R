# ==============================================================================
# Functions shared by the manuscript scripts
# "From post-war hardship to longevity leader: a period and cohort analysis of
#  the Spanish mortality paradox"
# ==============================================================================
library(tidyverse)
options(scipen = 999)

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
                                        ymin = 1930, ymax = 2019,
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
                                        ymin = 1940, ymax = 2019,
                                        age_starts = seq(0, 100, 10), width = 10,
                                        sexes = c("m", "f", "t"),
                                        mark_years = 1940:1960){

  cohs <- tibble()
  summ <- tibble()
  for(a in age_starts){
    for(ct in cts_ref){
      cat(paste0("ages ", a, "-", a + width - 1, " ", c1, " vs ", ct, "\n"))
      res <-
        est_coh_cont_diff_countries(c1 = c1, c2 = ct,
                                    amin = a, amax = a + width - 1,
                                    ymin = ymin, ymax = ymax,
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
# that crossed the age group during the marked period (per2).
# ------------------------------------------------------------------------------
plot_scc <- function(cohs,
                     cols = c("#AE2012", "#94D2BD", "#005F73", "#EE9B00"),
                     legend_title = "Gap between Spain and:",
                     drop_ags = "100_109"){

  cohs %>%
    filter(sex != "t", !(ags %in% drop_ags)) %>%
    mutate(sex = ifelse(sex == "f", "Female", "Male"),
           ags = str_replace(as.character(ags), "_", "-"),
           ags = fct_reorder(ags, amin),
           cts = str_remove(cts, "^[A-Z_]+?_")) %>%
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
