# ==============================================================================
# 12. Supplementary Figure S6. Lexis surfaces of the annual changes in the
#     age-specific contributions to the life expectancy gap between each
#     reference country and Spain, 1940-2019, ages 0-99, by sex.
#
# Each cell is DD_x(y, y+1) = D_{x,y+1} - D_{x,y}: the change, from year y to
# y + 1, in the (symmetric Arriaga) contribution of age x to the gap in life
# expectancy at birth between the reference country and Spain (REF - ESP).
# Negative values (blue) are changes in Spain's favour (a narrowing of the gap
# while Spain is behind), positive values (red) changes against Spain. The
# cohort contributions of Figure 5 are sums of these cells along Lexis
# diagonals, within 10-year age windows; here the whole surface is shown so
# that period blocks (vertical bands) and cohort traces (diagonals) can be
# told apart directly. Colours are assigned by quantiles of the negative and
# of the positive values pooled over all panels, so that the small
# contributions of adult ages remain visible next to the large infant ones.
#
# Input : data_inter/ltables_hmd.rds (01_hmd_life_tables.R)
# Output: data_inter/lexis_gap_contributions.rds
#         figures/supplement/figS6_lexis_gap_contributions.png
# ==============================================================================
rm(list = ls()); gc()
source("code/00_functions.R")

cts_ref <- c("SWE", "CHE", "NOR", "NLD")
names_cts <- c(SWE = "Sweden", CHE = "Switzerland", NOR = "Norway", NLD = "Netherlands")
ymin <- 1940
ymax <- 2019
amin <- 0
amax <- 99

hmd_lt <-
  read_rds("data_inter/ltables_hmd.rds") %>%
  filter(code %in% c("ESP", cts_ref))

# annual Arriaga contributions of each age to the gap REF - ESP ------------------
contr <- tibble()
for (ct in cts_ref) {
  cat(ct, "\n")
  for (s in c("f", "m")) {
    for (y in ymin:ymax) {
      tmp <- est_age_cont(c1 = "ESP", s1 = s, y1 = y,
                          c2 = ct, s2 = s, y2 = y,
                          amin = amin, amax = amax)
      contr <- bind_rows(contr, tmp %>% select(age, dt) %>%
                           mutate(cts = ct, sex = s, year = y))
    }
  }
}

# change from year y to y + 1 at each age -----------------------------------------
dd <-
  contr %>%
  arrange(cts, sex, age, year) %>%
  group_by(cts, sex, age) %>%
  mutate(ddt = lead(dt) - dt) %>%
  ungroup() %>%
  filter(!is.na(ddt)) %>%
  mutate(cohort = year - age,
         sex = ifelse(sex == "f", "Female", "Male"),
         country = factor(names_cts[cts], levels = names_cts))

write_rds(dd, "data_inter/lexis_gap_contributions.rds")

# colour bins: quantiles of the negative and the positive values (pooled) --------
probs <- c(0, .5, .8, .9, .95, .98, 1)
neg <- -dd$ddt[dd$ddt < 0]
pos <-  dd$ddt[dd$ddt > 0]
br_neg <- -rev(quantile(neg, probs))
br_pos <- quantile(pos, probs)
breaks <- unique(c(br_neg[1] * 1.001, br_neg[-1], br_pos[-1] * c(rep(1, length(br_pos) - 2), 1.001)))
# labels of the bins
lab <- sprintf("%.3f", breaks)
bin_labels <- paste0("(", lab[-length(lab)], ", ", lab[-1], "]")
cols <- c(rev(c("#deebf7", "#9ecae1", "#6baed6", "#3182bd", "#08519c", "#08306b")),
          c("#fee0d2", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d", "#99000d"))

dd2 <-
  dd %>%
  mutate(bin = cut(ddt, breaks = breaks, labels = bin_labels, include.lowest = TRUE))

stopifnot(!any(is.na(dd2$bin)))

# figure -------------------------------------------------------------------------
dd2 %>%
  ggplot(aes(year, age)) +
  facet_grid(sex ~ country) +
  geom_tile(aes(fill = bin)) +
  scale_fill_manual(values = setNames(cols, bin_labels), name = "Change in the age\ncontribution to the\ngap REF - ESP (years)",
                    drop = FALSE) +
  # Lexis grid and the 1940-1960 period
  geom_vline(xintercept = seq(1950, 2010, 10), linewidth = .15, linetype = "dashed", colour = "grey30") +
  geom_hline(yintercept = seq(10, 90, 10), linewidth = .15, linetype = "dashed", colour = "grey30") +
  geom_abline(intercept = seq(-2010, -1850, 10), slope = 1, linewidth = .15, linetype = "dashed", colour = "grey30") +
  annotate("rect", xmin = ymin - .5, xmax = 1960.5, ymin = amin - .5, ymax = amax + .5,
           fill = NA, colour = "black", linewidth = .4) +
  scale_x_continuous(breaks = seq(1940, 2020, 20), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 100, 20), expand = c(0, 0)) +
  coord_equal() +
  labs(x = "Period", y = "Age") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold", size = 9),
        panel.grid = element_blank(),
        panel.spacing = unit(.4, "cm"),
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.key.size = unit(.4, "cm"),
        plot.background = element_rect(fill = "white", colour = NA))

dir.create("figures/supplement", showWarnings = FALSE)
ggsave("figures/supplement/figS6_lexis_gap_contributions.png", w = 10, h = 5.6, dpi = 300)
