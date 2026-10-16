# ==============================================================================
# Exploratory plots of the cohort contributions to the life expectancy gap
# (relative, absolute and negative contributions by age group). Not in the
# manuscript. Requires data_inter/cohort_conts_gap_changes.RData produced by
# code/07_fig5_cohort_contributions_gap.R
# ==============================================================================
rm(list = ls()); gc()
source("code/exploratory/00_functions_exploratory.R")
ymin <- 1940; ymax <- 2019
dir.create("figures/exploratory/cohort_cont_gap", showWarnings = FALSE)

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

ggsave(paste0("figures/exploratory/cohort_cont_gap/cohort_rel_cont_gap_change_", ymin, "_", ymax, ".pdf"),
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

ggsave(paste0("figures/exploratory/cohort_cont_gap/cohort_abs_cont_gap_change_", ymin, "_", ymax, ".pdf"),
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

ggsave(paste0("figures/exploratory/cohort_cont_gap/cohort_rel_cont_gap_change_", ymin, "_", ymax, "2.png"),
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

ggsave(paste0("figures/exploratory/cohort_cont_gap/cohort_neg_cont_gap_change_", ymin, "_", ymax, "2.png"),
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

ggsave(paste0("figures/exploratory/cohort_cont_gap/cohort_coh_cont_stand_gap_change_", ymin, "_", ymax, "2.png"),
       w = 7, h = 7,
       limitsize = FALSE)


dt_gaps3 %>% 
  group_by(sex, ags, cts) %>% 
  mutate(dt_coh_sd = dt_coh/max(abs(dt_coh))) %>% 
  ungroup()
