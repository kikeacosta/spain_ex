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

y <- 1940
out1 <- tibble()
for(y in 1940:2019){
  
  tmp1 <- 
    est_age_cont(c1 = "ESP", s1 = "f", y1 = y-1, 
                 c2 = "ESP", s2 = "f", y2 = y) %>% 
    select(age, dd) %>% 
    mutate(year = y,
           code = "ESP")
  
  tmp2 <- 
    est_age_cont(c1 = "SWE", s1 = "f", y1 = y-1, 
                 c2 = "SWE", s2 = "f", y2 = y) %>% 
    select(age, dd) %>% 
    mutate(year = y,
           code = "SWE")
  
  out1 <- 
    out1 %>% 
    bind_rows(tmp1,
              tmp2)
  
}

# age-specific contributions to ex changes over time
out1 %>% 
  filter(age %in% 15:90,
         year %in% 1960:2019) %>% 
  ggplot()+
  facet_grid(~code)+
  geom_tile(aes(x = year, y = age, fill = dd))+
  scale_fill_gradient2(low = "red",
                       mid = "white",
                       high = "blue",
                       midpoint = 0)+
  coord_fixed()

dt2 <- 
  out1 %>% 
  group_by(code, age) %>% 
  mutate(dt_ch = dd - lag(dd))

# changes in age-specific contributions to ex changes over time
dt2 %>% 
  filter(age %in% 15:90,
         year %in% 1960:2019) %>% 
  ggplot()+
  facet_grid(~code)+
  geom_tile(aes(x = year, y = age, fill = dt_ch))+
  scale_fill_gradient2(low = "red",
                       mid = "white",
                       high = "blue",
                       midpoint = 0)+
  scale_y_continuous(breaks = seq(0, 100, 10))+
  coord_fixed()


# difference between Spain and Sweden in age-specific contributions to ex 
# changes over time
diffs <- 
  out1 %>% 
  spread(code, dd) %>% 
  mutate(ct_diff = ESP - SWE)


diffs %>% 
  filter(age %in% 15:90,
         year %in% 1960:1990) %>% 
  ggplot()+
  geom_tile(aes(x = year, y = age, fill = ct_diff))+
  scale_fill_gradient2(low = "red",
                       mid = "white",
                       high = "blue",
                       midpoint = 0)+
  coord_fixed()


# changes in the difference between Spain and Sweden in age-specific 
# contributions to ex changes over time
dt3 <- 
  diffs %>% 
  group_by(age) %>% 
  mutate(diff_ch = ct_diff - lag(ct_diff))

dt3 %>% 
  filter(age %in% 15:90,
         year %in% 1960:1990) %>% 
  ggplot()+
  geom_tile(aes(x = year, y = age, fill = diff_ch))+
  scale_fill_gradient2(low = "red",
                       mid = "white",
                       high = "blue",
                       midpoint = 0)+
  coord_fixed()

