rm (list = ls()); gc()
source("code/exploratory/00_functions_exploratory.R")

hmd_lt <- read_rds("data_inter/ltables_hmd.rds")

ct <- "ESP"
ymn <- 1970
ymx <- 2019
amn <- 50
amx <- 80

ctr <- "ESP"
# {
res <- 
  est_coh_cont_per(c = ctr, 
                   amin = 60, amax = 69,
                   ymin = 1960, ymax = 2019)

p1 <- 
  res$summ %>% 
  ggplot()+
  geom_point(aes(cohort, dt, col = sex))+
  geom_line(aes(cohort, dt, col = sex))+
  geom_hline(yintercept = 0, lty = "dashed")+
  scale_x_continuous(breaks = seq(1910, 1950, 5))+
  labs(y = "Total contribution\nto ex change (years)",
       x = "Cohort",
       color = "Sex",
       title = ctr)+
  theme_bw()+
  theme(strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

p2 <- 
  res$long %>% 
  group_by(cohort, sex) %>% 
  mutate(dt_av = mean(dt)) %>% 
  ggplot()+
  geom_point(aes(cohort, dt, col = sex), alpha = 0.1)+
  geom_violin(aes(cohort, dt, col = sex, group = interaction(sex, cohort)))+
  geom_point(aes(cohort, dt_av, col = sex), alpha = 1)+
  geom_hline(yintercept = 0, lty = "dashed")+
  scale_x_continuous(breaks = seq(1910, 1950, 5))+
  scale_color_discrete(labels = c("Females", "Males", "Total"))+
  # facet_grid(~sex)+
  labs(y = "Annual contributions\nto ex change (years)",
       x = "Cohort",
       color = "Sex")+
  theme_bw()+
  theme(strip.background = element_blank(),
        legend.position = "bottom")

# p1+p2
plot_grid(p1, p2, ncol = 1, rel_heights = c(.8, 1))

# ggsave(paste0("figures/ext_abs/coh_contrbs_abs_all_", ctr, ".png"),
#        w = 8, h = 5)
# }


res_esp <- est_coh_cont_per(c = "ESP", 
                            amin = 50, amax = 80,
                            ymin = 1970, ymax = 2019)

res_swe <- est_coh_cont_per(c = "SWE", 
                            amin = 50, amax = 80,
                            ymin = 1970, ymax = 2019)

res_che <- est_coh_cont_per(c = "CHE", 
                            amin = 50, amax = 80,
                            ymin = 1970, ymax = 2019)

res_esp$long

esp <- 
  res_esp$summ %>% 
  select(cohort, sex, dt_esp = dt)

swe <- 
  res_swe$summ %>% 
  select(cohort, sex, dt_swe = dt)

che <- 
  res_che$summ %>% 
  select(cohort, sex, dt_swe = dt)

dt <- 
  left_join(esp, swe, che) %>% 
  mutate(diff = dt_esp - dt_swe)

dt %>% 
  group_by(sex) %>% 
  summarise(dt_esp = sum(dt_esp),
            dt_swe = sum(dt_swe),
            diff = sum(diff))

dt %>% 
  ggplot()+
  geom_point(aes(cohort, diff, col = sex))+
  geom_line(aes(cohort, diff, col = sex))+
  geom_hline(yintercept = 0, lty = "dashed")+
  scale_x_continuous(breaks = seq(1910, 1950, 5))+
  labs(y = "Total contribution\nto ex change (years)",
       x = "Cohort",
       color = "Sex",
       title = "cohort contributions to ex change diff ESP and SWE")+
  theme_bw()+
  theme(strip.background = element_blank(),
        # axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")



dt <- 
  bind_rows(res_esp$summ,
            res_swe$summ,
            res_che$summ)

dt %>% 
  ggplot()+
  facet_grid(~sex)+
  geom_point(aes(cohort, dt, col = code))+
  geom_line(aes(cohort, dt, col = code))+
  geom_hline(yintercept = 0, lty = "dashed")+
  scale_x_continuous(breaks = seq(1910, 1950, 5))+
  labs(y = "Total contribution\nto ex change (years)",
       x = "Cohort",
       color = "Sex",
       title = "cohort contributions to ex change")+
  theme_bw()+
  theme(strip.background = element_blank(),
        # axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom"
        )


a1 <- 40
a2 <- 59
a <- 0
dt2 <- tibble()
for(a in seq(0, 100, 10)){
  swe_esp <- 
    est_coh_cont_diff_countries(c1 = "ESP", c2 = "SWE", 
                                amin = a, amax = a+9,
                                ymin = 1930, ymax = 2019)
  
  che_esp <- 
    est_coh_cont_diff_countries(c1 = "ESP", c2 = "CHE", 
                                amin = a, amax = a+9,
                                ymin = 1930, ymax = 2019)
  
  # dnk_esp <- 
  #   est_coh_cont_diff_countries(c1 = "ESP", c2 = "DNK", 
  #                               amin = a, amax = a+9,
  #                               ymin = 1930, ymax = 2019)
  
  nor_esp <- 
    est_coh_cont_diff_countries(c1 = "ESP", c2 = "NOR", 
                                amin = a, amax = a+9,
                                ymin = 1930, ymax = 2019)
  
  nld_esp <- 
    est_coh_cont_diff_countries(c1 = "ESP", c2 = "NLD", 
                                amin = a, amax = a+9,
                                ymin = 1930, ymax = 2019)
  
  dt2 <- 
    bind_rows(swe_esp$summ,
              che_esp$summ,
              nor_esp$summ,
              nld_esp$summ)
  
  # dt2 %>% 
  #   ggplot()+
  #   facet_grid(~sex)+
  #   geom_point(aes(cohort, dt, col = cts))+
  #   geom_line(aes(cohort, dt, col = cts))+
  #   geom_hline(yintercept = 0, lty = "dashed")+
  #   scale_x_continuous(breaks = seq(1800, 2000, 5))+
  #   labs(y = "Total contribution\nto ex Gap (years)",
  #        x = "Cohort",
  #        color = "Sex",
  #        title = paste0("cohort contributions to ex Gap, ages: ", a, "-", a+9))+
  #   theme_bw()+
  #   theme(strip.background = element_blank(),
  #         # axis.text.x = element_blank(),
  #         axis.title.x = element_blank(),
  #         legend.position = "bottom"
  #   )
  # ggsave(paste0("figures/cohort_cont_gap/cohort_cont_gap_ages_", a, "_", a+9, ".png"),
  #        w = 16,
  #        h = 8)
  
}

swe_esp <- 
  est_coh_cont_diff_countries(c1 = "ESP", c2 = "SWE", 
                              amin = a1, amax = a2,
                              ymin = 1930, ymax = 2019)

che_esp <- 
  est_coh_cont_diff_countries(c1 = "ESP", c2 = "CHE", 
                              amin = a1, amax = a2,
                              ymin = 1930, ymax = 2019)

dnk_esp <- 
  est_coh_cont_diff_countries(c1 = "ESP", c2 = "DNK", 
                              amin = a1, amax = a2,
                              ymin = 1930, ymax = 2019)

dt2 <- 
  bind_rows(swe_esp$summ,
            che_esp$summ,
            dnk_esp$summ)

dt2 %>% 
  ggplot()+
  facet_grid(~sex)+
  geom_point(aes(cohort, dt, col = cts))+
  geom_line(aes(cohort, dt, col = cts))+
  geom_hline(yintercept = 0, lty = "dashed")+
  scale_x_continuous(breaks = seq(1800, 2000, 5))+
  labs(y = "Total contribution\nto ex Gap (years)",
       x = "Cohort",
       color = "Sex",
       title = paste0("cohort contributions to ex Gap, ages: ", a1, "-", a2))+
  theme_bw()+
  theme(strip.background = element_blank(),
        # axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom"
  )
ggsave(paste0("figures/cohort_cont_gap/cohort_cont_gap_ages_", a1, "_", a2, ".png"),
       w = 16,
       h = 8)

