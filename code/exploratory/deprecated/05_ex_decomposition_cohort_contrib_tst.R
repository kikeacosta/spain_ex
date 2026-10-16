rm (list = ls()); gc()
source("code/exploratory/00_functions_exploratory.R")

hmd_lt <- read_rds("data_inter/ltables_hmd.rds")

hmd_lt %>% 
  filter(year %in% c(1940, 1960),
         sex == "t",
         age == 0,
         code %in% c("ESP", "SWE", "NOR", "CHE", "NLD")) %>% 
  select(code, year, ex) %>% 
  spread(year, ex)
  

dt2 <- tibble()
a <- 0

swe_esp <- 
  est_coh_cont_diff_countries(c1 = "ESP", c2 = "SWE", 
                              amin = a, amax = a+110,
                              ymin = 1940, ymax = 2019)


nor_esp <- 
  est_coh_cont_diff_countries_tst(c1 = "ESP", c2 = "NOR", 
                              amin = a, amax = a+110,
                              ymin = 1940, ymax = 2019)

tst <- 
  swe_esp %>% 
  group_by(sex, age) %>% 
  arrange(sex, year, age) %>% 
  mutate(dt_ch = dt - lag(dt)) %>% 
  ungroup()

copy_this(tst)

ages <- 
  tst %>% 
  drop_na(dt_ch) %>% 
  reframe(dt_ch = sum(dt_ch),
          .by = c(sex, age))

copy_this(ages %>% 
            filter(sex == "t"))

chts <- 
  tst %>% 
  drop_na(dt_ch) %>% 
  reframe(dt_ch = sum(dt_ch),
          ns = n(),
          .by = c(sex, cohort))

copy_this(chts %>% 
            filter(sex == "t"))




d_temp_sx <- 
  est_age_cont(c1 = "ESP", s1 = "t", y1 = 2000, 
               c2 = "SWE", s2 = "t", y2 = 2000,
               amin = amin, amax = amax) %>% 
  select(age, dd, di, dt) %>% 
  mutate(sex = s,
         year = i)





dt2 <- 
  tst %>% 
  reframe(dt_ch = sum(dt_ch),
          .by = c(sex, cohort))

copy_this(dt2)

dt2 <- 
  swe_esp$summ


est_coh_cont_diff_countries_tst <- 
  function(c1 = "ESP", c2 = "SWE", 
           amin = 10, amax = 85,
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
      mutate(cohort = year - age) %>%
      group_by(sex) %>% 
      mutate(yrs = paste0(ymin, "_", ymax),
             ags = paste0(amin, "_", amax),
             chs = paste0(min(cohort), "_", max(cohort)),
             cts = paste0(c1, "_", c2))
    
    return(dt1)
  }



swe_esp <- 
  est_coh_cont_diff_countries(c1 = "ESP", c2 = "SWE", 
                              amin = 0, amax = 110,
                              ymin = 1930, ymax = 2019)


swe_esp <- 
  est_coh_cont_diff_countries_tst(c1 = "ESP", c2 = "SWE", 
                                  amin = 0, amax = 110,
                                  ymin = 1930, ymax = 2019)

tst <- 
  swe_esp %>% 
  group_by(sex, age) %>% 
  arrange(sex, year, age) %>% 
  mutate(dt_ch = dt - lag(dt)) %>% 
  ungroup()

copy_this(tst)
