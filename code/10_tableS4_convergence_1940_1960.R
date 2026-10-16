# ==============================================================================
# 10. Supplementary Table S4 (and numbers quoted in the Results). Life
#     expectancy at birth in Spain and the gap with each reference country and
#     with their unweighted mean, in 1940, 1950 and 1960, by sex: relative
#     change of Spanish life expectancy and relative reduction of the gap over
#     the 1940s, the 1950s and 1940-1960.
#
# Input : data_inter/e0_hmd.rds (01_hmd_life_tables.R)
# Output: data_inter/tableS4_convergence_1940_1960.csv
# ==============================================================================
rm(list = ls()); gc()
library(tidyverse)

e0 <- read_rds("data_inter/e0_hmd.rds")

cts_ref <- c("NLD", "NOR", "SWE", "CHE")
yrs <- c(1940, 1950, 1960)

le <-
  e0 %>%
  filter(Code %in% c("ESP", cts_ref), Year %in% yrs) %>%
  pivot_longer(c(Female, Male, Total), names_to = "sex", values_to = "ex") %>%
  select(year = Year, code = Code, sex, ex)

# reference: each country and the unweighted mean of the four
ref <-
  le %>%
  filter(code != "ESP") %>%
  bind_rows(le %>% filter(code != "ESP") %>%
              reframe(ex = mean(ex), .by = c(year, sex)) %>%
              mutate(code = "MEAN")) %>%
  rename(ref = code, ex_ref = ex)

tab <-
  le %>%
  filter(code == "ESP") %>%
  select(year, sex, ex_esp = ex) %>%
  inner_join(ref, by = c("year", "sex")) %>%
  mutate(gap = ex_esp - ex_ref) %>%
  arrange(sex, ref, year) %>%
  group_by(sex, ref) %>%
  reframe(ex_esp_1940 = ex_esp[year == 1940],
          ex_esp_1950 = ex_esp[year == 1950],
          ex_esp_1960 = ex_esp[year == 1960],
          gap_1940 = gap[year == 1940],
          gap_1950 = gap[year == 1950],
          gap_1960 = gap[year == 1960],
          # relative change of Spanish life expectancy (%)
          esp_chg_1940s = 100 * (ex_esp_1950 / ex_esp_1940 - 1),
          esp_chg_1950s = 100 * (ex_esp_1960 / ex_esp_1950 - 1),
          esp_chg_1940_60 = 100 * (ex_esp_1960 / ex_esp_1940 - 1),
          # relative reduction of the gap (%), each decade relative to its start
          gap_red_1940s = 100 * (1 - gap_1950 / gap_1940),
          gap_red_1950s = 100 * (1 - gap_1960 / gap_1950),
          gap_red_1940_60 = 100 * (1 - gap_1960 / gap_1940)) %>%
  mutate(ref = factor(ref, levels = c(cts_ref, "MEAN")),
         sex = factor(sex, levels = c("Total", "Female", "Male"))) %>%
  arrange(sex, ref) %>%
  mutate(across(where(is.numeric), ~ round(.x, 1)))

print(as.data.frame(tab))

write_csv(tab, "data_inter/tableS4_convergence_1940_1960.csv")
