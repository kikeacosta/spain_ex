rm (list = ls()); gc()
source("code/exploratory/00_functions_exploratory.R")
library("tidyverse")

hmd_lt <- read_rds("data_inter/ltables_hmd.rds")

unique(hmd_lt$code)

lt1 <- 
  hmd_lt %>% 
  filter(code == "ESP",
         year == 1950,
         sex == "t")

lt2 <- 
  hmd_lt %>% 
  filter(code == "ESP",
         year == 2000,
         sex == "t")

c1 = "ESP"; s1 = "t"; y1 = 1938; 
c2 = "ESP"; s2 = "t"; y2 = 1939

esp_1910_1960 <- 
  est_age_cont(c1 = "ESP", s1 = "t", y1 = 1910, 
               c2 = "ESP", s2 = "t", y2 = 1960)

esp_1960_2010 <- 
  est_age_cont(c1 = "ESP", s1 = "t", y1 = 1960, 
               c2 = "ESP", s2 = "t", y2 = 2010)

esp_1910_1960 %>% summarise(dt = sum(dt))
esp_1910_1960 %>% filter(age == 0) %>% mutate(dt2 = ex2-ex1) %>% pull(dt2)

esp_1910_1960 %>% 
  ggplot()+
  geom_point(aes(age, dt))+
  theme_bw()+
  scale_x_continuous(breaks = seq(0, 110, 10))

esp_1960_2010 %>% 
  ggplot()+
  geom_point(aes(age, dt))+
  theme_bw()+
  scale_x_continuous(breaks = seq(0, 110, 10))

# i <- 1910

c_esp <- est_coh_cont_all(c = "ESP")

# selecting cohorts and ages to look at
c_esp2 <- 
  c_esp %>% 
  filter(cohort %in% 1910:1950,
         age %in% 20:69)

# changes by year
c_year <- 
  c_esp2 %>% 
  filter(age == 20) %>% 
  select(year, dt_all, dd_all, di_all)

# aggregating contributions from each cohort
c_esp3 <- 
  c_esp2 %>% 
  group_by(cohort) %>% 
  summarise(dd = sum(dd),
            di = sum(di),
            dt = sum(dt),
            )






c_ita <- est_coh_cont(c = "ITA")
c_che <- est_coh_cont(c = "CHE")

coh_com <- 
  bind_rows(c_esp,
            c_ita,
            c_che)

coh_com %>% 
  filter(cohort >= 1900) %>% 
  ggplot()+
  geom_point(aes(cohort, cont, col = code))+
  geom_line(aes(cohort, cont, col = code))+
  theme_bw()+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_x_continuous(breaks = seq(1900, 2020, 10))

coh_ratios <- 
  coh_com %>% 
  spread(code, cont) %>% 
  mutate(esp_ita = ESP-ITA,
         esp_che = ESP-CHE)

coh_ratios %>% 
  filter(cohort >= 1880) %>% 
  ggplot()+
  geom_point(aes(cohort, esp_ita), col = "red")+
  geom_point(aes(cohort, esp_che), col = "blue")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_bw()+
  scale_x_continuous(breaks = seq(1880, 2020, 10))


