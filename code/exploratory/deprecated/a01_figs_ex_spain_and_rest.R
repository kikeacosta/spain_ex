rm (list = ls()); gc()
# library("HMDHFDplus")
source("code/exploratory/00_functions_exploratory.R")
library("tidyverse")

hmd_lt <- read_rds("data_inter/ltables_hmd.rds")


unique(hmd_lt$code) %>% sort

# countries in 1920
cds <- 
  hmd_lt %>% 
  filter(year == 1920) %>% 
  pull(code) %>% 
  unique()

# annual life expectancy ranking
rnks <- 
  hmd_lt %>% 
  filter(age == 0) %>% 
  group_by(year, sex) %>% 
  arrange(-ex) %>% 
  mutate(rk = 1:n()) %>% 
  select(code, year, sex, ex, rk)

# annual life expectancy ranking among Western European countries
cts_wst <- c("ITA", "ESP", "FRATNP", "PRT", "GBR_NP", "BEL", "CHE",
             "NLD", "GRC", "AUT") %>% sort
cts_wst <- c("ITA", "ESP", "FRATNP", "PRT", "GBR_NP", "BEL", "CHE",
             "NLD", "GRC", "AUT", "AUS", "CAN", "USA", "NZL_NP",
             "NOR", "SWE", "DNK") %>% sort
rnks_wst <- 
  hmd_lt %>% 
  filter(age == 0,
         code %in% cts_wst) %>% 
  group_by(year, sex) %>% 
  arrange(-ex) %>% 
  mutate(rk = 1:n()) %>% 
  select(code, year, sex, ex, rk)

# Spain in ranking
rnks %>% 
  filter(code == "ESP") %>% 
  ggplot()+
  geom_point(aes(year, rk, col = sex), alpha = 0.6)+
  theme_bw()

rnks %>% 
  filter(code == "ESP") %>% 
  filter(year >= 1990) %>% 
  ggplot()+
  geom_point(aes(year, rk, col = sex), alpha = 0.6)+
  theme_bw()

# Italy vs Spain females
unique(rnks$code) %>% sort

rnks %>% 
  filter(code %in% c("ITA", "ESP", "FRATNP", "PRT", "GBR_NP", "BEL", "CHE",
                     "NLD", "GRC", "AUT"),
         sex == "f") %>% 
  filter(year %in% 1950:2010) %>% 
  ggplot()+
  geom_line(aes(year, rk, col = code), alpha = 1)+
  theme_bw()


# Western countries ====
# ~~~~~~~~~~~~~~~~~~~~~~
rnks %>% 
  filter(code %in% cts_wst,
         # sex == "f"
         ) %>% 
  filter(year %in% 1930:2019) %>% 
  ggplot()+
  facet_grid(~sex)+
  geom_line(aes(year, ex, col = code, linewidth = code), alpha = .7)+
  scale_x_continuous(breaks = seq(1930, 2020, 10))+
  # scale_linewidth_manual(values = c(rep(0.3, 3), .8, rep(0.3, 6)))+
  scale_linewidth_manual(values = c(rep(0.1, 6), .5, rep(0.1, 10)))+
  scale_alpha_manual(values = c(rep(0.2, 6), 1, rep(0.2, 10)))+
  guides(col = guide_legend(ncol = 10),
         linewidth = guide_legend(ncol = 10))+
  # labs(title = "Females")+
  theme_bw()+
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 6),
    legend.spacing = unit(0, 'cm'),
    legend.key.size = unit(0.4, 'cm'),
    axis.text.x = element_text(size = 5) 
  )


ggsave("figures/a_descriptive/ex_across_cts_wst_tst.png",
       w = 8,
       h = 4,
       dpi = 1000)

rnks %>% 
  filter(
    # code %in% c("ITA", "ESP", "FRATNP", "PRT", "GBR_NP", "BEL", "CHE",
    #                  "NLD", "GRC", "AUT"),
    # sex == "f"
    ) %>% 
  filter(year %in% 1930:2019) %>% 
  mutate(code = factor(code, levels = cts_ord)) %>% 
  ggplot()+
  facet_grid(~sex)+
  geom_line(aes(year, ex, col = code, linewidth = code), alpha = .7)+
  scale_x_continuous(breaks = seq(1930, 2020, 10))+
  scale_color_manual(values = c(rep("grey", 11), "black", rep("grey", 29)))+
  scale_linewidth_manual(values = c(rep(0.3, 11), .8, rep(0.3, 29)))+
  coord_cartesian(xlim = c(1930, 2020))+
  theme_bw()+
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 5) 
  )

ggsave("figures/a_descriptive/ex_across_cts_all.png",
       w = 8,
       h = 4,
       dpi = 600)


# all countries ====
# ~~~~~~~~~~~~~~~~~~
cts_ord <- unique(rnks$code) %>% sort

yrs <- c(1950:1985)
rnks %>% 
  filter(
    # sex == "f",
    year %in% yrs
  ) %>% 
  mutate(code = factor(code, levels = cts_ord)) %>% 
  ggplot()+
  facet_grid(~sex)+
  geom_line(aes(year, ex, col = code, linewidth = code), alpha = .7)+
  scale_x_continuous(breaks = seq(1930, 2020, 10))+
  scale_color_manual(values = c(rep("grey", 9), "black", rep("grey", 31)))+
  scale_linewidth_manual(values = c(rep(0.3, 9), .8, rep(0.3, 31)))+
  coord_cartesian(xlim = c(yrs[1], yrs[length(yrs)]))+
  # labs(title = "Females")+
  theme_bw()+
  theme(
    legend.text = element_text(size = 3),
    legend.title = element_text(size = 4),
    legend.spacing = unit(0, 'cm'),
    legend.key.size = unit(0.2, 'cm')
  )

ggsave(paste0("figures/a_descriptive/ex_across_cts_all_", 
              yrs[1], 
              "_", 
              yrs[length(yrs)],
              ".png"),
       w = 8,
       h = 4,
       dpi = 600)

# ranking ====
library(ggbump)

rnks_wst %>% 
  filter(year %in% c(seq(1930, 2015, 5), 2019)) %>% 
  ggplot(aes(year, rk, col = code))+
  facet_grid(~sex)+
  geom_point(alpha = 0.6)+
  geom_bump(aes(alpha = code), size = 1) +
  geom_text(data = rnks_wst %>% 
              filter(year == 1930),
            aes(x = year - 1, label = code),
            size = 2, hjust = 1) +
  geom_text(data = rnks_wst %>% 
              filter(year == 2019),
            aes(x = year + 1, label = code),
            size = 2, hjust = 0) +
  scale_y_continuous(breaks = 0:20)+
  scale_x_continuous(breaks = seq(1930, 2020, 10))+
  scale_alpha_manual(values = c(rep(0.3, 6), 1, rep(0.3, 10)))+
  coord_cartesian(xlim = c(1925, 2025))+
  labs(title = "rank in ex among Western countries")+
  theme_bw()+
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 8) 
  )

ggsave(paste0("figures/a_descriptive/ex_rank_wst.png"),
       w = 10, h = 5,
       dpi = 1000)


# changes in ex

ex_ch <- 
  rnks %>% 
  group_by(code, sex) %>% 
  arrange(year) %>% 
  mutate(ex_ch = ex/lag(ex))

cts_60_19 <- 
  ex_ch %>% 
  filter(year == 1960 | year == 2019, 
         sex == "t") %>% 
  group_by(code) %>% 
  filter(n() == 2) %>% 
  pull(code) %>% unique()

chunk <- 
  ex_ch %>% 
  filter(code == "ESP",
         sex == "t")

library(mgcv)

rnks_smt <- 
  rnks %>% 
  filter(year %in% 1919:2019) %>% 
  group_by(code, sex) %>% 
  do(smt_ex(chunk = .data)) %>% 
  do(smt_ex_rg(chunk = .data)) %>% 
  ungroup()

rnks_smt %>% 
  filter(code %in% cts_wst,
         year >= 1920) %>% 
  ggplot()+
  facet_grid(~sex)+
  geom_point(aes(year, ex, col = code), alpha = 0.1, size = .15)+
  geom_line(aes(year, ex_smt, col = code, linewidth = code, alpha = code))+
  scale_x_continuous(breaks = seq(1920, 2020, 10))+
  scale_y_continuous(breaks = seq(40, 90, 5))+
  scale_alpha_manual(values = c(rep(0.3, 6), 1, rep(0.3, 10)))+
  scale_linewidth_manual(values = c(rep(0.2, 6), .5, rep(0.2, 10)))+
  # labs(title = paste0("Sex_", sx))+
  guides(col = guide_legend(ncol = 10),
         linewidth = guide_legend(ncol = 10),
         alpha = guide_legend(ncol = 10))+
  theme_bw()+
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 6),
    legend.spacing = unit(0, 'cm'),
    legend.key.size = unit(0.4, 'cm'),
    axis.text.x = element_text(size = 5) 
  )

ggsave(paste0("figures/a_descriptive/ex_across_cts_wst_smt.png"),
       w = 8,
       h = 4,
       dpi = 1000)

rnks_smt %>% 
  filter(code %in% cts_wst,
         year >= 1920) %>% 
  ggplot()+
  facet_grid(~sex)+
  geom_point(aes(year, ex, col = code), alpha = 0.2, size = .3)+
  geom_line(aes(year, ex_smt_rg, col = code, linewidth = code, alpha = code))+
  scale_x_continuous(breaks = seq(1920, 2020, 10))+
  scale_y_continuous(breaks = seq(40, 90, 5))+
  scale_alpha_manual(values = c(rep(0.3, 6), 1, rep(0.3, 10)))+
  scale_linewidth_manual(values = c(rep(0.2, 6), .5, rep(0.2, 10)))+
  # labs(title = paste0("Sex_", sx))+
  guides(col = guide_legend(ncol = 10),
         linewidth = guide_legend(ncol = 10),
         alpha = guide_legend(ncol = 10))+
  theme_bw()+
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 6),
    legend.spacing = unit(0, 'cm'),
    legend.key.size = unit(0.4, 'cm'),
    axis.text.x = element_text(size = 5) 
  )

ggsave(paste0("figures/a_descriptive/ex_across_cts_wst_smt_rg.png"),
       w = 8,
       h = 4,
       dpi = 1000)


ex_ch <- 
  rnks_smt %>% 
  group_by(code, sex) %>% 
  arrange(year) %>% 
  mutate(ex_ch = ex/lag(ex),
         ex_smt_ch = ex_smt/lag(ex_smt),
         ex_smt_rg_ch = ex_smt_rg/lag(ex_smt_rg),
  )

ex_ch %>% 
  filter(code %in% cts_wst,
         # code != "PRT",
         # sex != "t",
         year %in% 1930:2019,
         # code %in% cts_60_19
  ) %>% 
  drop_na(ex_smt_ch) %>% 
  ggplot()+
  facet_grid(~sex)+
  geom_tile(aes(year, code, fill = ex_smt_ch))+
  scale_fill_gradient2(midpoint = 1)+
  scale_x_continuous(breaks = seq(1920, 2020, 10))+
  coord_cartesian(expand = 0)+
  theme(
    axis.text.x = element_text(size = 3) 
  )+
  theme_bw()

ggsave(paste0("figures/a_descriptive/heat_ex_smt_change_cts_wst.png"),
       w = 10,
       h = 5)

ex_ch %>% 
  filter(code %in% cts_wst,
         # code != "PRT",
         # sex != "t",
         year %in% 1930:2019,
         # code %in% cts_60_19
  ) %>% 
  drop_na(ex_smt_rg_ch) %>% 
  ggplot()+
  facet_grid(~sex)+
  geom_tile(aes(year, code, fill = ex_smt_rg_ch))+
  scale_fill_gradient2(midpoint = 1)+
  scale_x_continuous(breaks = seq(1920, 2020, 10))+
  coord_cartesian(expand = 0)+
  theme(
    axis.text.x = element_text(size = 3) 
  )+
  theme_bw()

ggsave(paste0("figures/a_descriptive/heat_ex_smt_change_cts_wst_rg.png"),
       w = 10,
       h = 5)





ex_ch %>% 
  filter(code %in% cts_wst,
         # code != "PRT",
         # sex != "t",
         year %in% 1960:2000,
         # code %in% cts_60_19
  ) %>% 
  drop_na(ex_smt_ch) %>% 
  ggplot()+
  facet_grid(~sex)+
  geom_tile(aes(year, code, fill = ex_smt_rg_ch))+
  scale_fill_gradient2(midpoint = 1)+
  scale_x_continuous(breaks = seq(1920, 2020, 10))+
  coord_cartesian(expand = 0)+
  theme(
    axis.text.x = element_text(size = 5) 
  )+
  theme_bw()

ggsave("figures/a_descriptive/heat_ex_smt_change_cts_wst_rg_60_00.png",
       w = 10,
       h = 5)


ex_ch %>% 
  filter(code %in% cts_wst,
         # sex != "t",
         # year >= 1920,
         year %in% 1930:2019) %>% 
  ggplot()+
  facet_grid(~sex)+
  # geom_point(aes(year, ex, col = code), alpha = 0.2, size = .8)+
  geom_line(aes(year, ex_smt_ch, col = code, linewidth = code,
                alpha = code))+
  scale_x_continuous(breaks = seq(1930, 2020, 10))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  # scale_y_continuous(breaks = seq(60, 90, 5))+
  scale_alpha_manual(values = c(rep(0.5, 6), 1, rep(0.5, 10)))+
  scale_linewidth_manual(values = c(rep(0.3, 6), .7, rep(0.3, 10)))+
  # labs(title = paste0("Sex_", sx))+
  guides(col = guide_legend(ncol = 10),
         linewidth = guide_legend(ncol = 10),
         alpha = guide_legend(ncol = 10))+
  theme_bw()+
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 6),
    legend.spacing = unit(0, 'cm'),
    legend.key.size = unit(0.4, 'cm'),
    axis.text.x = element_text(size = 5) 
  )

ggsave(paste0("figures/a_descriptive/ex_smt_change_cts_wst.png"),
       w = 10, h = 5,
       dpi = 600)



ch_rk_wst <- 
  ex_ch %>% 
  filter(code %in% cts_wst,
         year %in% 1930:2019) %>% 
  group_by(year, sex) %>% 
  arrange(-ex_smt_ch) %>% 
  mutate(rk = 1:n()) %>% 
  arrange(-ex_smt_rg_ch) %>% 
  mutate(rk_rg = 1:n()) 

rnks_smt
ch_rk_wst %>% 
  filter(year %in% c(seq(1930, 2015, 5), 2019)) %>% 
  ggplot(aes(year, rk, col = code))+
  facet_grid(~sex)+
  geom_point(alpha = 0.6)+
  geom_bump(aes(alpha = code), size = 1) +
  geom_text(data = ch_rk_wst %>% 
              filter(year == 1930),
            aes(x = year - 1, label = code),
            size = 2, hjust = 1) +
  geom_text(data = ch_rk_wst %>% 
              filter(year == 2019),
            aes(x = year + 1, label = code),
            size = 2, hjust = 0) +
  scale_y_continuous(breaks = 0:20)+
  scale_x_continuous(breaks = seq(1930, 2020, 10))+
  scale_alpha_manual(values = c(rep(0.3, 6), 1, rep(0.3, 10)))+
  coord_cartesian(xlim = c(1923, 2027))+
  labs(title = "rank in ex change among Western countries")+
  theme_bw()+
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 8) 
  )

ggsave(paste0("figures/a_descriptive/ex_change_rank_wst.png"),
       w = 10, h = 5,
       dpi = 1000)

ch_rk_wst %>% 
  filter(year %in% c(seq(1930, 2015, 5), 2019)) %>% 
  ggplot(aes(year, rk_rg, col = code))+
  facet_grid(~sex)+
  geom_point(alpha = 0.6)+
  geom_bump(aes(alpha = code), size = 1) +
  geom_text(data = ch_rk_wst %>% 
              filter(year == 1930),
            aes(x = year - 1, label = code),
            size = 2, hjust = 1) +
  geom_text(data = ch_rk_wst %>% 
              filter(year == 2019),
            aes(x = year + 1, label = code),
            size = 2, hjust = 0) +
  scale_y_continuous(breaks = 0:20)+
  scale_x_continuous(breaks = seq(1930, 2020, 10))+
  scale_alpha_manual(values = c(rep(0.3, 6), 1, rep(0.3, 10)))+
  coord_cartesian(xlim = c(1923, 2027))+
  labs(title = "rank in ex change among Western countries (super smoothed)")+
  theme_bw()+
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 8) 
  )

ggsave(paste0("figures/a_descriptive/ex_change_rigid_rank_wst.png"),
       w = 10, h = 5,
       dpi = 1000)



  