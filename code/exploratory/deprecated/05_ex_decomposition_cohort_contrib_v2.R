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



amin = 10; amax = 80
c1 <- "ESP"
y1 <- 2000
s1 <- "t"

ltt <- 
  hmd_lt %>% 
  filter(code == c1,
         year == y1,
         age %in% amin:amax,
         sex == s1)

T1 <- 
  ltt %>% 
  reframe(tx = sum(Lx)) %>% 
  pull(tx)
  
l1 <- 
  ltt %>% 
  filter(age == amin) %>% 
  pull(lx)

nex <- T1/l1

tx1 <- 
  ltt %>% 
  filter(age == amin) %>% 
  pull(Tx)

tx2 <- 
  ltt %>% 
  filter(age == amax) %>% 
  pull(Tx)

# re-estimate partial ex from complete life table
tx <- 
  ltt %>% 
  mutate(Lx = ifelse(age != max(age), 
                     lead(lx) * (lead(age) - age) + dx * ax, 
                     lx / mx),
         Tx = order_by(-age, cumsum(Lx)),
         ex = round(Tx/lx, 2))

Tx <- 
  tx %>% 
  filter(age == amin) %>% 
  pull(Tx)

lx <- 
  ltt %>% 
  filter(age == amin) %>% 
  pull(lx)

par_ex <- Tx/lx

c1 = "ESP"; s1 = "t"; y1 = 1938; amin = 10; amax = 85
c2 = "ESP"; s2 = "t"; y2 = 1970; 

lt1 <- 
  hmd_lt %>% 
  filter(code == c1,
         year == y1,
         age %in% amin:amax,
         sex == s1) %>% 
  mutate(Lx = ifelse(age != max(age), 
                     lead(lx) * (lead(age) - age) + dx * ax, 
                     lx / mx),
         Tx = order_by(-age, cumsum(Lx)),
         ex = round(Tx/lx, 2))

lt2 <- 
  hmd_lt %>% 
  filter(code == c2,
         year == y2,
         age %in% amin:amax,
         sex == s2) %>% 
  mutate(Lx = ifelse(age != max(age), 
                     lead(lx) * (lead(age) - age) + dx * ax, 
                     lx / mx),
         Tx = order_by(-age, cumsum(Lx)),
         ex = round(Tx/lx, 2))


e_01 <- 
  lt1 %>% 
  filter(age == amin) %>% 
  pull(ex)

e_02 <- 
  lt2 %>% 
  filter(age == amin) %>% 
  pull(ex)

e_diff <- e_02 - e_01
e_diff

l01 <-   
  lt1 %>% 
  filter(age == amin) %>% 
  pull(lx)


lt1_arr <- 
  lt1 %>% 
  select(age, lx1 = lx, Lx1 = Lx, Tx1 = Tx, ex1 = ex) %>% 
  left_join(lt2 %>% 
              select(age, lx2 = lx, Lx2 = Lx, Tx2 = Tx, ex2 = ex))

lt1_arr2 <- 
  lt1_arr %>% 
  mutate(dd = lx1/l01*(Lx2/lx2-Lx1/lx1),
         di = lead(Tx2)/l01*(lx1/lx2-lead(lx1)/lead(lx2)),
         dt = ifelse(age < amax, dd + di, dd))

lt1_arr2 %>% 
  summarise(dt = sum(dt))


