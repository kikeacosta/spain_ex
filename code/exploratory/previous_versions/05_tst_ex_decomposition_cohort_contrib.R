rm (list = ls()); gc()
source("code/00_functions.R")

hmd_lt <- read_rds("data_inter/ltables_hmd.rds")

ymin = 1940; ymax = 2019

dt_gaps_cohs <- tibble()
dt_gaps_long <- tibble()
dt_gaps_summ <- tibble()

for(a in seq(0, 100, 10)){
  swe_esp <- 
    est_coh_cont_diff_countries(c1 = "ESP", c2 = "SWE", 
                                amin = a, amax = a+9,
                                ymin = ymin, ymax = ymax)
  
  che_esp <- 
    est_coh_cont_diff_countries(c1 = "ESP", c2 = "CHE", 
                                amin = a, amax = a+9,
                                ymin = ymin, ymax = ymax)
  
  nor_esp <- 
    est_coh_cont_diff_countries(c1 = "ESP", c2 = "NOR", 
                                amin = a, amax = a+9,
                                ymin = ymin, ymax = ymax)
  
  nld_esp <- 
    est_coh_cont_diff_countries(c1 = "ESP", c2 = "NLD", 
                                amin = a, amax = a+9,
                                ymin = ymin, ymax = ymax)
  
  dt_gaps_summ <- 
    dt_gaps_summ %>% 
    bind_rows(swe_esp$summ,
              che_esp$summ,
              nor_esp$summ,
              nld_esp$summ)
  
  dt_gaps_cohs <- 
    dt_gaps_cohs %>% 
    bind_rows(swe_esp$coh_cont,
              che_esp$coh_cont,
              nor_esp$coh_cont,
              nld_esp$coh_cont)
  
  dt_gaps_long <- 
    dt_gaps_long %>% 
    bind_rows(swe_esp$long,
              che_esp$long,
              nor_esp$long,
              nld_esp$long)
  
}

dt_gaps2 <- 
  dt_gaps_cohs %>% 
  mutate(ags = factor(ags, levels = c("0_9", "10_19", "20_29", "30_39", "40_49", 
                                      "50_59", "60_69", "70_79", "80_89", 
                                      "90_99", "100_109")),
         age = ifelse(ags == "100_109", "100", str_sub(ags, 1, 2)),
         age = str_replace(age, "_", ""),
         age = as.integer(age) + 5,
         per = cohort + age,
         per_rap = ifelse(per %in% 1940:1960, 1, NA),
         per2 = cohort*per_rap)


save(dt_gaps2, dt_gaps_summ, file = "data_inter/cohort_conts_gap_changes.RData")


load("data_inter/cohort_conts_gap_changes.RData")

dt_gaps2 %>% 
  group_by(ags, sex, cts) %>% 
  summarise(dt_rel = sum(dt_rel))

dt_gaps2 %>% 
  ggplot()+
  facet_grid(ags~sex, scales = "free_y", space = "free_y")+
  geom_line(aes(cohort, dt_rel, col = cts), alpha = .7)+
  geom_hline(yintercept = 0, lty = "dashed")+
  scale_x_continuous(breaks = seq(1800, 2020, 10))+
  labs(y = "Total contribution\nto ex Gap reduction (%)",
       x = "Cohort",
       color = "Sex",
       title = paste0("cohort contributions to ex Gap reduction, period:", ymin, "-", ymax))+
  theme_bw()+
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "bottom"
  )

ggsave(paste0("figures/cohort_cont_gap/cohort_rel_cont_gap_change_", ymin, "_", ymax, ".pdf"),
       w = 16, h = 20,
       limitsize = FALSE)


dt_gaps2 %>% 
  ggplot()+
  facet_grid(ags~sex, scales = "free_y", space = "free_y")+
  geom_line(aes(cohort, dt_coh, col = cts), alpha = .7)+
  geom_hline(yintercept = 0, lty = "dashed")+
  scale_x_continuous(breaks = seq(1800, 2020, 10))+
  labs(y = "Total contribution\nto ex Gap reduction (years)",
       x = "Cohort",
       color = "Sex",
       title = paste0("cohort contributions to ex Gap reduction, period:", ymin, "-", ymax))+
  theme_bw()+
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "bottom"
  )

ggsave(paste0("figures/cohort_cont_gap/cohort_abs_cont_gap_change_", ymin, "_", ymax, ".pdf"),
       w = 16, h = 20,
       limitsize = FALSE)


dt_gaps3 <- 
  dt_gaps2 %>% 
  mutate(age = ifelse(ags == "100_109", "100", str_sub(ags, 1, 2)),
         age = str_replace(age, "_", ""),
         age = as.integer(age) + 5,
         per = cohort + age,
         per_rap = ifelse(per %in% 1940:1960, 1, NA),
         per2 = cohort*per_rap) %>% 
  group_by(sex, ags, cts) %>% 
  mutate(dt_coh_sd = dt_coh/max(abs(dt_coh))) %>% 
  ungroup()


unique(dt_gaps3$ags)
unique(dt_gaps3$age)

cols4 <- c("#AE2012", "#94D2BD", "#005F73", "#EE9B00")

dt_gaps3 %>% 
  filter(sex != "t") %>% 
  mutate(sex = ifelse(sex == "f", "Female", "Male")) %>% 
  ggplot()+
  facet_grid(sex~ags, scales = "free", space = "free")+
  geom_line(aes(cohort, dt_rel, col = cts), alpha = .7)+
  geom_point(aes(x = per2, y = Inf), col = "black")+
  geom_hline(yintercept = 0, lty = "dashed")+
  scale_x_continuous(breaks = seq(1800, 2020, 10))+
  scale_color_manual(values = cols4)+
  labs(y = "Total contribution to partial ex Gap reduction (%)",
       x = "Cohort",
       color = "Gap between:",
       # title = paste0("cohort contributions to ex Gap reduction, period:", ymin, "-", ymax)
       )+
  coord_flip()+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 7),
        # axis.title.x = element_blank(),
        legend.position = "bottom"
  )

ggsave(paste0("figures/manuscript/cohort_rel_cont_gap_change_", ymin, "_", ymax, "2.png"),
       w = 7, h = 7,
       limitsize = FALSE)

dt_gaps3 %>% 
  filter(sex != "t") %>% 
  mutate(sex = ifelse(sex == "f", "Female", "Male")) %>% 
  ggplot()+
  facet_grid(sex~ags, scales = "free", space = "free")+
  geom_line(aes(cohort, dt_neg, col = cts), alpha = .7)+
  geom_point(aes(x = per2, y = Inf), col = "black")+
  geom_hline(yintercept = 0, lty = "dashed")+
  scale_x_continuous(breaks = seq(1800, 2020, 10))+
  scale_color_manual(values = cols4)+
  labs(y = "Total contribution to partial ex Gap reduction (%)",
       x = "Cohort",
       color = "Gap between:",
       # title = paste0("cohort contributions to ex Gap reduction, period:", ymin, "-", ymax)
  )+
  coord_flip()+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 7),
        # axis.title.x = element_blank(),
        legend.position = "bottom"
  )

ggsave(paste0("figures/manuscript/cohort_neg_cont_gap_change_", ymin, "_", ymax, "2.png"),
       w = 7, h = 7,
       limitsize = FALSE)

dt_gaps3 %>% 
  filter(sex != "t") %>% 
  mutate(sex = ifelse(sex == "f", "Female", "Male")) %>% 
  ggplot()+
  facet_grid(sex~ags, scales = "free", space = "free")+
  geom_line(aes(cohort, dt_coh_sd, col = cts), alpha = .4, linewidth = 0.1)+
  geom_point(aes(cohort, dt_coh_sd, col = cts), alpha = .7, size = 0.1)+
  geom_point(aes(x = per2, y = Inf), col = "black")+
  geom_hline(yintercept = 0, lty = "dashed")+
  scale_x_continuous(breaks = seq(1800, 2020, 10))+
  scale_color_manual(values = cols4)+
  labs(y = "Standardized Contribution to Change in Partial \nLife Expectancy Gap (Max = 1 per Age Group)",
       x = "Cohort",
       color = "Gap between:",
       # title = paste0("cohort contributions to ex Gap reduction, period:", ymin, "-", ymax)
  )+
  coord_flip()+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 7),
        # axis.title.x = element_blank(),
        legend.position = "bottom"
  )

ggsave(paste0("figures/manuscript/cohort_coh_cont_stand_gap_change_", ymin, "_", ymax, "2.png"),
       w = 7, h = 7,
       limitsize = FALSE)


dt_gaps3 %>% 
  group_by(sex, ags, cts) %>% 
  mutate(dt_coh_sd = dt_coh/max(abs(dt_coh))) %>% 
  ungroup()
tst <- 
  dt_gaps3 %>% 
  select(cts, cohort, sex, ags, dt_coh_sd)

dt_gaps3 %>% 
  filter(sex != "t", ags != "100_109") %>% 
  mutate(sex = ifelse(sex == "f", "Female", "Male"),
         ags = str_replace(ags, "_", "-"),
         cts = str_remove(cts, "ESP_")) %>% 
  ggplot(aes(cohort, dt_coh_sd))+
  facet_grid(fct_rev(ags)~sex, scales = "free", space = "free", switch = "y")+
  geom_line(aes(cohort, dt_coh_sd, col = cts), alpha = .3, lwd = 0.1)+
  geom_point(aes(cohort, dt_coh_sd, col = cts), alpha = .5, size = 0.03)+
  geom_point(aes(x = per2, y = Inf), col = "#033b86", size = .5)+
  geom_smooth(aes(group = 1), method = "loess",  
              span = 0.4,          # DECREASED value for more flexibility/local fit
              se = F,
              lwd = 0.3, alpha = .6,
              color = "black") +
  geom_hline(yintercept = 0, lty = "dashed", lwd = .2)+
  scale_x_continuous(breaks = seq(1800, 2020, 10))+
  scale_y_continuous(position = "right",
                     breaks = c(-1, 0, 1))+
  scale_color_manual(values = cols4)+
  labs(y = "Standardized Contribution to LE gap",
       x = "Cohort",
       color = "Gap between Spain and:",
       # title = paste0("cohort contributions to ex Gap reduction, period:", ymin, "-", ymax)
  )+
  coord_cartesian(ylim = c(-1,1),
                  expand = 0)+
  # coord_flip()+
  guides(
    color = guide_legend(
      # Increase size of points and thickness of lines in the legend symbols
      override.aes = list(size = 2.5, linewidth = 1.5) 
    )
  )+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9),
        # strip.position = "left",
        # strip.text.y.left = element_text(angle = 0),
        panel.spacing.x = unit(0.5, "cm"),
        # panel.spacing.y = unit(0.3, "cm"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        # axis.text.y.right = element_text(size = 7),
        axis.title.y.right = element_text(vjust = 0.5, angle = 90),# axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(.8, "cm")
  )

ggsave(paste0("figures/manuscript/cohort_coh_cont_stand_gap_change_", ymin, "_", ymax, "tst.png"),
       w = 6, h = 6,
       dpi = 600,
       limitsize = FALSE)


