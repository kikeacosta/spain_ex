rm (list = ls()); gc()
source("code/exploratory/00_functions_exploratory.R")

hmd_lt <- read_rds("data_inter/ltables_hmd.rds")
dt <- read_rds("data_inter/deaths_exposures_hmd.rds")


unique(hmd_lt$code)

war <- 
  hmd_lt %>% 
  filter(code == "ESP",
         year %in% 1934:1943,
         sex == "t")

pwr <- 
  hmd_lt %>% 
  filter(code == "ESP",
         year %in% 1940:1959,
         sex == "t")

war %>% 
  ggplot()+
  geom_line(aes(age, mx, col = factor(year)))+
  scale_y_log10()+
  scale_x_continuous(breaks = seq(0, 100, 10))+
  theme_bw()

pwr %>% 
  ggplot()+
  geom_line(aes(age, mx, col = factor(year)))+
  scale_y_log10()+
  scale_x_continuous(breaks = seq(0, 100, 10))+
  theme_bw()

dt2 <- 
  dt %>% 
  filter(code == "ESP",
         year %in% 1920:2019,
         age %in% 10:100)


sxs <- unique(dt2$sex)
out <- tibble()

for(s in sxs){
  
  tmp <- 
    dt %>% 
    filter(code == "ESP",
           year %in% 1920:2019,
           age %in% 10:100,
           sex == s) %>% 
    mutate(w = ifelse(year %in% 1935:1960, 0, 1))
  
  # We need to give the information in matrix forms (age x period)
  ylist <- unique(tmp$year) %>% sort() # all periods in data, ordered
  alist <- unique(tmp$age) %>% sort() # all ages in data, ordered
  
  # mortality matrix
  deaths <- 
    matrix(tmp$dts, 
           nrow = length(alist), 
           ncol = length(ylist), 
           byrow = F)
  colnames(deaths) <- ylist
  rownames(deaths) <- alist
  
  # lets have a look of the deaths matrix 
  deaths
  
  # population at risk (population/exposure) matrix
  exposure <- 
    matrix(tmp$exposure, 
           nrow = length(alist),
           ncol = length(ylist), 
           byrow = F)
  colnames(exposure) <- ylist
  rownames(exposure) <- alist
  
  # weights for interpolation matrix
  ws <- 
    matrix(tmp$w, 
           nrow = length(alist),
           ncol = length(ylist), 
           byrow = F)
  colnames(ws) <- ylist
  rownames(ws) <- alist
  
  
  # smoothing mortality in 2D using p-splines
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  fit_2d <- 
    Mort2Dsmooth(x = alist, 
                 y = ylist, 
                 Z = deaths, 
                 W = ws,
                 offset = log(exposure),
                 overdispersion = TRUE, 
                 method = 1)
  
  # extracting the smoothed log rates 
  fit_2d$logmortality
  # transforming them to smoothed deaths
  mx_smooth <- 
    (exp(fit_2d$logmortality))
  
  # from matrix to tidy (long) form
  smt <- 
    mx_smooth %>% 
    # keeping the row names as age values
    as_tibble(rownames = "age") %>%
    mutate(age = age %>% as.double()) %>% 
    # reshaping to long form (tidy)
    gather(-age, key = year, value = mx) %>% 
    mutate(type = "m_smoothed",
           year = as.integer(year)) %>% 
    # replacing missing values and estimating log_rates (/100k)
    replace_na(list(mx = 0)) %>% 
    mutate(sex = s)

  out <- 
    out %>% 
    bind_rows(smt)
  
}


exc <- 
  dt2 %>% 
  select(year, sex, age, dts, exposure) %>% 
  left_join(out %>% 
              select(year, sex, age, mx_smt = mx)) %>% 
  filter(year %in% 1935:1960) %>% 
  mutate(bsn = mx_smt * exposure,
         cohort = year - age) %>% 
  reframe(dts = sum(dts),
          bsn = sum(bsn),
          .by = c(sex, cohort)) %>% 
  mutate(exc = dts - bsn,
         psc = dts / bsn,
  ) %>% 
  filter(cohort %in% 1910:1950)


exc %>% 
  filter(sex != "total") %>% 
  ggplot()+
  geom_point(aes(cohort, psc, col = sex))+
  geom_line(aes(cohort, psc, col = sex))+
  theme_bw()+
  scale_y_log10()+
  scale_x_continuous(breaks = seq(1900, 1960, 5))+
  geom_hline(yintercept = 1, lty = "dashed")+
  theme(legend.position = "bottom")

ggsave("figures/ext_abs/relative_excess_1935_1960.png",
       w = 8, h = 4)  
