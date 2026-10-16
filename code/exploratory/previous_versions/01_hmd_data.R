rm (list = ls()); gc()
library("HMDHFDplus")
library("tidyverse")

# obtaining all country codes from the HMD 
cds_hmd <- getHMDcountries()[3] %>% pull()
# cds_hmd2 <- cds_hmd[!cds_hmd %in% cts_exclude]

# HMD user and password
# ~~~~~~~~~~~~~~~~~~~~
hmd_us <- Sys.getenv("hmd_us")
hmd_pw <- Sys.getenv("hmd_pw")

# identifying those with data for 2020
hmd <- tibble()
for(ct in cds_hmd){
  
  cat(paste0(ct, "\n"))
  chunk_d <- 
    readHMDweb(ct, "E0per", hmd_us, hmd_pw) %>%
    # filter(Year >= 2010) %>%
    as_tibble() %>%
    mutate(Code = ct)
  
  hmd <- 
    hmd %>%
    bind_rows(chunk_d)
}


write_rds(hmd, "data_inter/e0_hmd.rds")
unique(hmd$Code)

cts_1920 <- 
  hmd %>% 
  group_by(Code) %>% 
  filter(min(Year) <= 1920) %>% 
  pull(Code) %>% 
  unique()

hmd %>% 
  filter(
    # Code %in% c("ESP", "ITA", "GBR_NP", "FRATNP", "GRC", "USA", "AUS"),
    Year >= 1900
  ) %>% 
  ggplot()+
  geom_line(aes(Year, Total, col = Code))+
  theme_bw()

hmd %>% 
  filter(
    Code %in% cts_1920,
    Year %in% 1920:2019
  ) %>% 
  ggplot()+
  geom_line(aes(Year, Total, col = Code))+
  theme_bw()

hmd %>% 
  filter(
    Code %in% "ESP",
    Year %in% 1920:2019
  ) %>% 
  ggplot()+
  geom_line(aes(Year, Total, col = Code))+
  theme_bw()


# loading complete life tables from HMD
hmd_lt <- tibble()
for(ct in cds_hmd){
  
  cat(paste0(ct, "\n"))
  
  # total sex
  chunk_b <- 
    readHMDweb(ct, "bltper_1x1", hmd_us, hmd_pw) %>%
    filter(Year >= 1900) %>%
    as_tibble() %>%
    mutate(code = ct,
           sex = "t") %>% 
    rename(year = Year,
           age = Age)
  
  # females
  chunk_f <- 
    readHMDweb(ct, "fltper_1x1", hmd_us, hmd_pw) %>%
    filter(Year >= 1900) %>%
    as_tibble() %>%
    mutate(code = ct,
           sex = "f") %>% 
    rename(year = Year,
           age = Age)

  # males
  chunk_m <- 
    readHMDweb(ct, "mltper_1x1", hmd_us, hmd_pw) %>%
    filter(Year >= 1900) %>%
    as_tibble() %>%
    mutate(code = ct,
           sex = "m") %>% 
    rename(year = Year,
           age = Age)
  
  hmd_lt <- 
    hmd_lt %>%
    bind_rows(chunk_b,
              chunk_f,
              chunk_m)
}

write_rds(hmd_lt, "data_inter/ltables_hmd.rds")

