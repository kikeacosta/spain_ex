# ==============================================================================
# Functions used only by the exploratory scripts (not needed for the manuscript)
# Sources the manuscript functions first.
# ==============================================================================
source("code/00_functions.R")
library(MortalitySmooth)
library(cowplot)
library(mgcv)

copy_this <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,file = paste0("clipboard-", object.size(x)),sep="\t",row.names=row.names,col.names=col.names,...)
}

# Using the Arriaga decomposition to obtain cohort-specific contributions to 
# life expectancy changes
est_coh_cont_per <- function(c = "ESP", 
                         amin = 10, amax = 85,
                         ymin = 1930, ymax = 2019){
  
  diffs <- tibble()
  for(i in ymin:(ymax-1)){
    
    for (s in c("m", "f", "t")){
      
      d_temp_sx <- 
        est_age_cont(c1 = c, s1 = s, y1 = i, 
                     c2 = c, s2 = s, y2 = i+1,
                     amin = amin, amax = amax) %>% 
        select(age, dd, di, dt) %>% 
        mutate(sex = s,
               year = i + 1)
      
      diffs <- 
        bind_rows(diffs, d_temp_sx)
      
    }
  }
  
  dt1 <- 
    diffs %>% 
    mutate(cohort = year - age) %>%
    group_by(cohort, sex) %>% 
    mutate(yrs = n()) %>% 
    ungroup() %>% 
    filter(yrs == max(yrs))
  
  dt2 <- 
    dt1 %>% 
    group_by(sex) %>% 
    mutate(pos = ifelse(dt > 0 , dt, 0), # looking at positive contributions
           neg = ifelse(dt < 0 , dt, 0), # looking at negative contributions
           pos_tot = sum(pos), # adding together positives
           neg_tot = sum(neg), # adding together positives
           dt_rel = pos/pos_tot) %>% # calculating the relative contribution to positive changes
    ungroup() %>% 
    mutate(yrs = paste0(ymin, "_", ymax),
           ags = paste0(amin, "_", amax),
           chs = paste0(min(cohort), "_", max(cohort)))
  
  coh_cont <- 
    dt2 %>% 
    reframe(dd = sum(dd),
            di = sum(di),
            dt = sum(dt),
            dt_pos = sum(pos),
            dt_neg = sum(neg),
            dt_rel = sum(dt_rel),
            .by = c(cohort, sex)) %>% 
    arrange(sex, cohort) %>% 
    mutate(yrs = paste0(ymin, "_", ymax),
           ags = paste0(amin, "_", amax),
           chs = paste0(min(cohort), "_", max(cohort)),
           code = c)
  
  out <- list(long = dt2, summ = coh_cont)
  return(out)
}

est_age_cont_diff_countries <- function(c1 = "ESP", c2 = "SWE", 
                                        amin = 0, amax = 110,
                                        ymin = 1930, ymax = 2019){
  
  diffs <- tibble()
  for(i in ymin:ymax){
    
    for (s in c("m", "f", "t")){
      
      d_temp_sx <- 
        est_age_cont(c1 = c1, s1 = s, y1 = i, 
                     c2 = c2, s2 = s, y2 = i,
                     amin = amin, amax = amax) %>% 
        select(age, dd, di, dt) %>% 
        mutate(sex = s,
               year = i)
      
      diffs <- 
        bind_rows(diffs, d_temp_sx)
      
    }
  }
  
  dt1 <- 
    diffs %>% 
    reframe(dd = sum(dd),
            di = sum(di),
            dt = sum(dt),
            .by = c(sex, age)) %>% 
    arrange(sex, age) %>% 
    mutate(yrs = paste0(ymin, "_", ymax),
           ags = paste0(amin, "_", amax),
           cts = paste0(c1, "_", c2))
  
  return(dt1)
}


# lexis surfaces of mortality change
plot_change <- function(c, s, amin, amax, ymin, ymax){
  # filtering data, and computing rates and log_rates 
  db2 <- 
    hmd %>% 
    mutate(deaths = dts + 1,
           Mx = 100000 * dts / exposure, 
           log_m = log(Mx)) %>% 
    filter(code == c,
           sex == s,
           age %in% amin:amax, 
           year %in% ymin:ymax)
  
  amin2 <- db2 %>% pull(age) %>% min()
  amax2 <- db2 %>% pull(age) %>% max()
  ymin2 <- db2 %>% pull(year) %>% min()
  ymax2 <- db2 %>% pull(year) %>% max()
  
  ylist <- unique(db2$year) %>% sort() # all periods in data, ordered
  alist <- unique(db2$age) %>% sort() # all ages in data, ordered
  
  # mortality
  deaths <- 
    matrix(db2$dts, nrow = length(alist), ncol = length(ylist), byrow = F)
  colnames(deaths) <- ylist
  rownames(deaths) <- alist
  
  # population at risk (population/exposure)
  exposure <- 
    matrix(db2$exposure, nrow = length(alist), ncol = length(ylist), byrow=F)
  colnames(exposure) <- ylist
  rownames(exposure) <- alist
  
  # smoothing mortality with best AIC (that is method=2, also possible best BIC with method=1)
  fit <- 
    Mort2Dsmooth(x = alist, y = ylist, Z = deaths, offset = log(exposure),
                 overdispersion = TRUE, method = 2)
  
  # transforming them to smoothed deaths
  mx_smooth <- (exp(fit$logmortality) * 100000)
  
  # from matrix to tidy form (adjusting ages to start in 0)
  smt <- 
    mx_smooth %>% 
    as_tibble(rownames = "age") %>%
    gather(-age, key = year, value = Mx) %>% 
    mutate(type = "m_smoothed",
           age = as.integer(age),
           year = as.integer(year))
  
  # replacing missing values and estimating log_rates (/100k)
  smt2 <- 
    smt %>% 
    replace_na(list(Mx = 0)) %>% 
    mutate(log_m = log(Mx))
  
  db_per <- 
    smt2 %>% 
    group_by(age) %>% 
    mutate(ch = ((Mx / lag(Mx)) - 1) * 100) %>% 
    ungroup() %>% 
    drop_na()
  
  qt <- 25
  # one for decrease of mortality (from green to blue), 
  # another for mortality increases (from yellow to red)
  col_scale <- c(colorRampPalette(c("royalblue", "springgreen"), space = "Lab")(qt),
                 colorRampPalette(c("yellow", "red"), space = "Lab")(qt))
  
  # Definition of brackets for the scale of change
  val <- unique(db_per$ch)
  # separate negative and positive values
  pval <- val[val>0] # all positive values of change (mortality deterioration)
  nval <- val[val<0] # all negative values of change (mortality improvement)
  # identification of brackets for positive values (minimum + 23 quantiles + maximum)
  pcop <-c(min(pval), quantile(pval, prob=1/(qt-1)*(1:(qt-2))), max(pval)*1.01)
  # the same as above but for negative values (minimum + 23 quantiles + maximum)
  ncop <-c(min(nval)*1.01, quantile(nval, prob=1/(qt-1)*(1:(qt-2))), max(nval)*1.01) 
  # chain of brackets 25 ranges for negative changes, central value of no change (0), 
  # and 25 ranges for positive changes
  breaks_mc <- c(ncop, 0, pcop) 
  
  # adding to each value of change (continuous) the corresponding bracket (a discrete interval)
  db_per2 <- 
    db_per %>% 
    mutate(ch_cut = cut(db_per$ch, breaks = breaks_mc)) 
  
  # all categories
  cuts <- db_per2 %>% pull(ch_cut) %>% unique() %>% sort()
  
  # assigning color to each category
  col_values <- setNames(col_scale, cuts)
  
  # Plot of mortality change over periods
  p_change <- 
    db_per2 %>% 
    ggplot(aes(year, age, z = ch_cut)) +
    geom_tile(aes(fill = ch_cut))+
    # adding the color palette constructed above
    scale_fill_manual(values = col_values, 
                      breaks = cuts[c(1, 10, 20, 25, 30, 40, 50)], 
                      name = "Mortality\nchange %")+ 
    #adding contour lines when the slope is 0
    geom_contour(aes(z = ch), breaks = 0, col="black", 
                 alpha=.7, linewidth = .3)+ 
    scale_x_continuous(expand = c(0,0), breaks = seq(ymin2, ymax2, 10)) +
    scale_y_continuous(expand = c(0,0), breaks = seq(amin2, amax2, 10))+
    # adding the grid of the Lexis diagram
    geom_vline(xintercept = seq(ymin2, ymax2, 10), linetype = "dashed", 
               color = "grey30", linewidth = .20, alpha = 0.8) +
    geom_hline(yintercept = seq(amin2, amax2, 10), linetype = "dashed", 
               color = "grey30", linewidth = .20, alpha = 0.8) +
    geom_rect(xmin = 2019.65, xmax = 2032, ymin = -1, ymax = 100.5, fill="white")+
    geom_rect(xmin = 1919, xmax = 2030, ymin = 100.05, ymax = 110.5, fill="white")+
    geom_abline(intercept = seq(-2020, -(ymin2-amax2), 10), slope = 1, 
                linetype = "dashed", color = "grey30", 
                linewidth = .2, alpha = 0.8)+
    labs(x="Period", y="Age")+
    # aesthetic details
    coord_equal() +
    labs(title = paste0(c, "_", s),
         x="Period", y="Age")+
    theme_minimal()+
    # deleting default discrete legend
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      plot.title = element_text(size = 12),
      panel.background = element_rect(fill = NA),
      panel.grid.major = element_line(colour = NA),
      panel.grid.minor = element_line(colour = NA),
      plot.background = element_rect(fill = "white", colour = "transparent")
    )
  
  return(p_change)
}




# lexis surfaces of mortality change
mort_change_data <- function(c, s, amin, amax, ymin, ymax){
  # filtering data, and computing rates and log_rates 
  db2 <- 
    hmd %>% 
    mutate(deaths = dts + 1,
           Mx = 100000 * dts / exposure, 
           log_m = log(Mx)) %>% 
    filter(code == c,
           sex == s,
           age %in% amin:amax, 
           year %in% ymin:ymax)
  
  amin2 <- db2 %>% pull(age) %>% min()
  amax2 <- db2 %>% pull(age) %>% max()
  ymin2 <- db2 %>% pull(year) %>% min()
  ymax2 <- db2 %>% pull(year) %>% max()
  
  ylist <- unique(db2$year) %>% sort() # all periods in data, ordered
  alist <- unique(db2$age) %>% sort() # all ages in data, ordered
  
  # mortality
  deaths <- 
    matrix(db2$dts, nrow = length(alist), ncol = length(ylist), byrow = F)
  colnames(deaths) <- ylist
  rownames(deaths) <- alist
  
  # population at risk (population/exposure)
  exposure <- 
    matrix(db2$exposure, nrow = length(alist), ncol = length(ylist), byrow=F)
  colnames(exposure) <- ylist
  rownames(exposure) <- alist
  
  # smoothing mortality with best AIC (that is method=2, also possible best BIC with method=1)
  fit <- 
    Mort2Dsmooth(x = alist, y = ylist, Z = deaths, offset = log(exposure),
                 overdispersion = TRUE, method = 2)
  
  # transforming them to smoothed deaths
  mx_smooth <- (exp(fit$logmortality) * 100000)
  
  # from matrix to tidy form (adjusting ages to start in 0)
  smt <- 
    mx_smooth %>% 
    as_tibble(rownames = "age") %>%
    gather(-age, key = year, value = Mx) %>% 
    mutate(type = "m_smoothed",
           age = as.integer(age),
           year = as.integer(year))
  
  # replacing missing values and estimating log_rates (/100k)
  smt2 <- 
    smt %>% 
    replace_na(list(Mx = 0)) %>% 
    mutate(log_m = log(Mx))
  
  db_per <- 
    smt2 %>% 
    group_by(age) %>% 
    mutate(ch = ((Mx / lag(Mx)) - 1) * 100) %>% 
    ungroup() %>% 
    drop_na()
  
  qt <- 25
  # one for decrease of mortality (from green to blue), 
  # another for mortality increases (from yellow to red)
  col_scale <- c(colorRampPalette(c("royalblue", "springgreen"), space = "Lab")(qt),
                 colorRampPalette(c("yellow", "red"), space = "Lab")(qt))
  
  # Definition of brackets for the scale of change
  val <- unique(db_per$ch)
  # separate negative and positive values
  pval <- val[val>0] # all positive values of change (mortality deterioration)
  nval <- val[val<0] # all negative values of change (mortality improvement)
  # identification of brackets for positive values (minimum + 23 quantiles + maximum)
  pcop <-c(min(pval), quantile(pval, prob=1/(qt-1)*(1:(qt-2))), max(pval)*1.01)
  # the same as above but for negative values (minimum + 23 quantiles + maximum)
  ncop <-c(min(nval)*1.01, quantile(nval, prob=1/(qt-1)*(1:(qt-2))), max(nval)*1.01) 
  # chain of brackets 25 ranges for negative changes, central value of no change (0), 
  # and 25 ranges for positive changes
  breaks_mc <- c(ncop, 0, pcop) 
  
  # adding to each value of change (continuous) the corresponding bracket (a discrete interval)
  db_per2 <- 
    db_per %>% 
    mutate(ch_cut = cut(db_per$ch, breaks = breaks_mc)) 
}

plot_mort_change_data <- function(dt){

  amin2 <- dt %>% pull(age) %>% min()
  amax2 <- dt %>% pull(age) %>% max()
  ymin2 <- dt %>% pull(year) %>% min()
  ymax2 <- dt %>% pull(year) %>% max()
  
  qt <- 25
  # one for decrease of mortality (from green to blue),
  # another for mortality increases (from yellow to red)
  col_scale <- c(colorRampPalette(c("royalblue", "springgreen"), space = "Lab")(qt),
                 colorRampPalette(c("yellow", "red"), space = "Lab")(qt))

  # all categories
  cuts <- dt %>% pull(ch_cut) %>% unique() %>% sort()

  # assigning color to each category
  col_values <- setNames(col_scale, cuts)

  # Plot of mortality change over periods
  p_change <-
    dt %>%
    ggplot(aes(year, age, z = ch_cut)) +
    geom_tile(aes(fill = ch_cut))+
    # adding the color palette constructed above
    scale_fill_manual(values = col_values,
                      breaks = cuts[c(1, 10, 20, 25, 30, 40, 50)],
                      name = "Mortality\nchange %")+
    #adding contour lines when the slope is 0
    geom_contour(aes(z = ch), breaks = 0, col="black",
                 alpha=.7, linewidth = .3)+
    scale_x_continuous(expand = c(0,0), breaks = seq(ymin2, ymax2, 10)) +
    scale_y_continuous(expand = c(0,0), breaks = seq(amin2, amax2, 10))+
    # adding the grid of the Lexis diagram
    geom_vline(xintercept = seq(ymin2, ymax2, 10), linetype = "dashed",
               color = "grey30", linewidth = .20, alpha = 0.8) +
    geom_hline(yintercept = seq(amin2, amax2, 10), linetype = "dashed",
               color = "grey30", linewidth = .20, alpha = 0.8) +
    geom_rect(xmin = 2019.65, xmax = 2032, ymin = -1, ymax = 100.5, fill="white")+
    geom_rect(xmin = 1919, xmax = 2030, ymin = 100.05, ymax = 110.5, fill="white")+
    geom_abline(intercept = seq(-2020, -(ymin2-amax2), 10), slope = 1,
                linetype = "dashed", color = "grey30",
                linewidth = .2, alpha = 0.8)+
    labs(x="Period", y="Age")+
    # aesthetic details
    coord_equal() +
    labs(title = paste0(c, "_", s),
         x="Period", y="Age")+
    theme_minimal()+
    # deleting default discrete legend
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      plot.title = element_text(size = 12),
      panel.background = element_rect(fill = NA),
      panel.grid.major = element_line(colour = NA),
      panel.grid.minor = element_line(colour = NA),
      plot.background = element_rect(fill = "white", colour = "transparent")
    )

  return(p_change)

}



get_drift <- 
  function(mdl){
    
    pcoeff <- 
      coef(summary(mdl))[,1] %>%
      as_tibble(rownames = "coeff") %>%
      filter(coeff %in% c("P", "C")) %>%
      pull(value)
    
    # this is in log scale, let's transform it in percentage value
    # change of mortality between periods (%)
    drift <- (exp(pcoeff) - 1)*100
    drift
    
    # interpretation
    cat(paste0(round(drift, 3),
               "% change in mortality each period/cohort category"))
  }

# extracting coefficients from an APC model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
extract_coeffs <- function(mod){
  tp1 <- 
    coef(summary(mod)) %>% 
    as_tibble(rownames = "coeff") %>% 
    mutate(tdim = case_when(str_detect(coeff, "\\(A\\)") ~ "Age",
                            str_detect(coeff, "\\(P\\)") ~ "Period",
                            str_detect(coeff, "\\(C\\)") ~ "Cohort",
                            str_detect(coeff, "I\\(C") ~ "Drift",
                            str_detect(coeff, "I\\(P") ~ "Drift"),
           effect = exp(Estimate),
           ll = exp(Estimate - 1.96*`Std. Error`),
           ul = exp(Estimate + 1.96*`Std. Error`)) %>% 
    separate(coeff, c("trash", "value"), sep = "\\)") %>% 
    mutate(value = value %>% as.double()) %>% 
    select(tdim, value, effect, ll, ul)
  
  ps_model <- tp1 %>% filter(tdim == "Period") %>% pull(value)
  ps_data <- unique(dt2$P)
  ps_miss <- ps_data[!(ps_data %in% ps_model)]
  
  cs_model <- tp1 %>% filter(tdim == "Cohort") %>% pull(value)
  cs_data <- unique(dt2$C)
  cs_miss <- cs_data[!(cs_data %in% cs_model)]
  
  tp2 <- 
    tp1 %>% 
    bind_rows(tibble(tdim = "Period", value = ps_miss, effect = 1, ll = 1, ul = 1),
              tibble(tdim = "Cohort", value = cs_miss, effect = 1, ll = 1, ul = 1)) %>% 
    mutate(tdim = factor(tdim, levels = c("Age", "Period", "Cohort"))) %>% 
    arrange(tdim, value) %>% 
    filter(tdim != "Drift")
  
  drift <- 
    tp1 %>% 
    filter(tdim == "Drift") %>% pull(effect)
  
  out <- 
    list("coeffs" = tp2, "drift" = drift)
  return(out)
}

# plotting Carstensen APC model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_carst <- 
  function(mod){
    
    age_eff <- 
      mod$Age %>% 
      as_tibble() %>% 
      mutate(dim = "Age")
    
    per_eff <- 
      mod$Per %>% 
      as_tibble() %>% 
      rename(Year = 1,
             RR = 2) %>% 
      mutate(dim = "Period")
    
    coh_eff <- 
      mod$Coh %>% 
      as_tibble() %>% 
      rename(Year = 1,
             RR = 2) %>% 
      mutate(dim = "Cohort")
    
    bind_rows(coh_eff, per_eff)
    
    p_a_eff <- 
      age_eff %>% 
      ggplot()+
      geom_line(aes(Age, Rate))+
      scale_y_log10()+
      facet_wrap(~dim, scales = "free_x")+
      theme_bw()
    
    p_pc_eff <- 
      bind_rows(coh_eff, per_eff) %>% 
      mutate(dim = factor(dim, levels = c("Period", "Cohort"))) %>% 
      ggplot()+
      geom_ribbon(aes(Year, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.3)+
      geom_line(aes(Year, RR))+
      scale_y_log10(breaks = c(0.2, 0.5, 0.7, 0.8, 1, 1.2, 1.5, 2, 4, 5))+
      scale_x_continuous(breaks = seq(1800, 2020, 10))+
      theme_bw()+
      facet_grid(~dim, scales = "free_x", space = "free_x")+
      geom_hline(yintercept = 1, linetype = "dashed")
    
    p_a_eff+p_pc_eff+ 
      plot_layout(widths = c(1, 2.5))
  }

smt_ex <- function(chunk){
  
  model <- 
    gam(ex ~ 
          s(year, bs = 'ps', m = c(2,2), k = 25), 
        # weights = w,
        data = chunk, 
        family = poisson(link = "log")
        , gamma = 0.01
    )
  
  res <- 
    predict(model, 
            newdata = chunk,
            type = "response", 
            se.fit = TRUE)
  
  chunk2 <- 
    chunk %>% 
    mutate(ex_smt = res$fit,
           ex_smt_l = res$fit - 1.96 * res$se.fit,
           ex_smt_u = res$fit + 1.96 * res$se.fit)
  
}

smt_ex_rg <- function(chunk){
  
  model <- 
    gam(ex ~ 
          s(year, bs = 'ps', m = c(2,2), k = 25), 
        # weights = w,
        data = chunk, 
        family = poisson(link = "log")
        , gamma = 0.5
    )
  
  res <- 
    predict(model, 
            newdata = chunk,
            type = "response", 
            se.fit = TRUE)
  
  chunk2 <- 
    chunk %>% 
    mutate(ex_smt_rg = res$fit,
           ex_smt_rg_l = res$fit - 1.96 * res$se.fit,
           ex_smt_rg_u = res$fit + 1.96 * res$se.fit)
  
}

