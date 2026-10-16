rm(list=ls())
source("code/exploratory/00_functions_exploratory.R")
library(Epi)
library(patchwork)
library(tidyverse)

hmd <- read_rds("data_inter/deaths_exposures_hmd.rds")

# let's have a look at the mortality experience of the US males since the 1950s
cd <- "ESP"
sx <- "total"

amin <- 10
amax <- 85
pmin <- 1960
pmax <- 2019

dt <- 
  hmd %>%
  filter(code == cd,
         sex == sx,
         age %in% amin:amax,
         year %in% pmin:pmax)

unique(dt$age)
unique(dt$year)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Descriptive APC analysis ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# to simplify the analysis we reduce the categories by
# grouping ages and periods in 5-year groups 
# (cohorts (P-A) will be automatically grouped too)
dt2 <- 
  dt %>% 
  rename(A = age,
         P = year) %>% 
  # creating new age categories grouped in 5 years 
  mutate(A = A - A%%5,
         P = P - P%%5,
         # creating a variable for birth cohort
         C = P - A) %>% 
  # aggregating deaths and exposures by age and period intervals
  summarise(dts = sum(dts) %>% round(),
            exposure = sum(exposure) %>% round(),
            .by = c(A, P, C)) %>% 
  # estimating death rates (/100K)
  mutate(mx = 1e5*dts/exposure) %>% 
  # removing weird rates
  filter(!is.na(mx) & mx!= Inf & mx!= 0)

dt2

# 4 classical plots for descriptive APC 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# In APC literature it is often suggested to build exploratory plots of the 
# linking the three dimensions Age-Period-Cohort

# a. Rates versus age at death for different periods:
# rates in the same age-group connected.
p1 <- 
  dt2 %>% 
  ggplot()+
  geom_line(aes(A, mx, col = factor(P), group = P))+
  scale_y_log10()+
  theme_bw()
p1
# b. Rates versus age at death for different birth cohorts:
# rates in the same birth-cohort connected.
p2 <- 
  dt2 %>% 
  ggplot()+
  geom_line(aes(A, mx, col = factor(C), group = C))+
  scale_y_log10()+
  theme_bw()
p2

p1+p2
# c. Rates versus date of death:
# rates in the same age-group connected by period.
p3 <- 
  dt2 %>% 
  ggplot()+
  geom_line(aes(P, mx, col = factor(A), group = A))+
  scale_y_log10()+
  theme_bw()
p3
# d. Rates versus date of date of birth:
# rates in the same age-group connected by cohort.
p4 <- 
  dt2 %>% 
  ggplot()+
  geom_line(aes(C, mx, col = factor(A), group = A))+
  scale_y_log10()+
  theme_bw()
p4

(p1+p2)/(p3+p4)

# ggsave("figures/basic_apc_desc_plots.png",
#        w = 10, h = 10)
# ~~~~~~~~~~~~~~~~~~~~~~


# now, stat APC models...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Statistical APC analysis ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# In the Epi package we need to rename the variables as 
# A: age
# P: year
# D: deaths
# Y: exposures

dt_carst <- 
  dt2 %>% 
  rename(D = dts,
         Y = exposure) %>% 
  select(-mx, -C)

# fitting the model
# ~~~~~~~~~~~~~~~~~
apc_c <- 
  apc.fit(dt_carst, 
          model = "factor", 
          dr.extr = "1", 
          parm = "AdPC", 
          scale = 10^5)

# model, effects, drift, APC fitting evaluation, etc.
apc_c

# drift
# ~~~~~
apc_c$Drift

# plotting
# ~~~~~~~~
# apc.plot(apc_c)
# I don't like the Epi function for plotting APC estimates
# then, I made my own function for plotting the APC estimates from Carstensen
# models
plot_carst(apc_c)

# now longitudinal age-specific death rates 
acp_c <- 
  apc.fit(dt_carst, 
          model = "factor", 
          dr.extr = "1", 
          parm = "AdCP", 
          scale = 10^5)

plot_carst(acp_c)

# the function is very versatile and has many ways to parameterize
# for instance, not grouping ages and periods, 

# age in single-year of period and age
dt_carst_x1 <- 
  dt %>% 
  # renaming variables for the Epi package
  select(D = dts,
         Y = exposure,
         A = age,
         P = year)

# for instance, we can add all the drift to the cohort dimension
acp_factor <- 
  apc.fit(dt_carst_x1, 
          model = "factor", 
          dr.extr = "Y", 
          # all drift to the C
          parm = "AdCP", 
          scale = 10^5)


a <- plot_carst(acp_factor)

# or we can fit splines for the nonlinear effects instead of using categorical 
# variables. We obtain smoothed APC effects, which are much more convenient 
# for analyses of changes in the trend
acp_splines <- 
  apc.fit(dt_carst_x1, 
          model = "bs", 
          # ref.c = 1950,
          # defining the amount of knots in each APC dimension for fitting the 
          # splines 
          # more knots means more flexibility
          npar = c(A = 10, P = 10, C = 15),
          dr.extr = "1", 
          parm = "AdCP", 
          scale = 10^5)

b <- plot_carst(acp_splines)
a/b
ggsave("figures/apc_effects_spain_total_ages_50_100.png",
       w = 14, h = 10)
