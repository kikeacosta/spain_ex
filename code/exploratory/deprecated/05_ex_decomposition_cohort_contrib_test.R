rm (list = ls()); gc()
source("code/exploratory/00_functions_exploratory.R")

hmd_lt <- read_rds("data_inter/ltables_hmd.rds")

# c1 = "ESP"; s1 = "m"; y1 = 1960; 
# c2 = "ESP"; s2 = "m"; y2 = 1961;
# amin = 0; amax = 110

sx <- "f"
ct <- "ESP"
ymn <- 1960
ymx <- 1980
amn <- 0
amx <- 60

dt <- 
  est_coh_cont(c = ct, s = sx,
               amin = amn, amax = amx,
               ymin = ymn, ymax = ymx)

# testing consistency
dt %>%
  summarise(sum(dt))

hmd_lt %>%
  filter(code == ct,
         sex == sx,
         age >= amn, 
         age <= amx,
         year == ymn | year == ymx) %>% 
  group_by(year) %>% 
  mutate(Tx = order_by(-age, cumsum(Lx)),
         ex = ifelse(lx>0, round(Tx/lx, 2), 0)) %>% 
  ungroup() %>% 
  filter(age == amn) %>% 
  mutate(per = ifelse(year == min(year), "per1", "per2")) %>% 
  select(per, ex) %>%
  spread(per, ex) %>%
  mutate(diff = per2 - per1)

coh_cmp <- 
  dt %>% 
  mutate(cohort = year - age) %>%
  group_by(cohort) %>% 
  mutate(yrs = n()) %>% 
  ungroup() %>% 
  filter(yrs == max(yrs))

coh_tot <- 
  coh_cmp %>% 
  reframe(dd = sum(dd),
          di = sum(di),
          dt = sum(dt),
          dd_rel = sum(dd_rel),
          di_rel = sum(di_rel),
          dt_rel = sum(dt_rel),
          .by = cohort)

coh_cnt_abs <- 
  coh_tot %>% 
  select(-dd_rel, -di_rel, -dt_rel) %>% 
  gather(-cohort, key = type, value = cont)

coh_cnt_rel <- 
  coh_tot %>% 
  select(-dd, -di, -dt) %>% 
  gather(-cohort, key = type, value = cont)

coh_cnt_rel2 <- 
  coh_tot %>% 
  select(-dd_rel, -di_rel, -dt_rel) %>% 
  gather(-cohort, key = type, value = cont) %>% 
  mutate(cont = ifelse(cont > 0, cont, 0)) %>% 
  group_by(cohort, type) %>% 
  summarise(con = sum(cont)) %>% 
  group_by(type) %>% 
  mutate(con_rel = 100*con/sum(con))


coh_cnt_abs %>% 
  filter(cohort %in% 1910:1939) %>% 
  ggplot()+
  geom_point(aes(cohort, cont, col = type))+
  geom_line(aes(cohort, cont, col = type))+
  geom_hline(yintercept = 0, lty = "dashed")+
  scale_x_continuous(breaks = seq(1910, 1950, 5))+
  scale_color_discrete(labels = c("Direct", "Indirect", "Total"))+
  labs(y = "Cont to ex change (years)",
       x = "Cohort",
       color = "Effect")+
  theme_bw()

ggsave(paste0("figures/ext_abs/coh_contrbs_abs_", sx, ".png"),
       w = 6, h = 3)

coh_cnt_rel2 %>% 
  filter(cohort %in% 1910:1939) %>% 
  ggplot()+
  geom_point(aes(cohort, con_rel, col = type))+
  geom_line(aes(cohort, con_rel, col = type))+
  scale_color_discrete(labels = c("Direct", "Indirect", "Total"))+
  scale_x_continuous(breaks = seq(1910, 1950, 5))+
  labs(y = "Cont to positive ex changes (%)",
       x = "Cohort",
       color = "Effect")+
  theme_bw()

ggsave(paste0("figures/ext_abs/coh_contrbs_rel_", sx, ".png"),
       w = 6, h = 3)
