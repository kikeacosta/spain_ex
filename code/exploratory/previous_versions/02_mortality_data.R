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
ct <- "AUS"
for(ct in cds_hmd){
  
  cat(paste0(ct, "\n"))
  chunk_d <- 
    readHMDweb(ct, "Deaths_1x1", hmd_us, hmd_pw) %>%
    # filter(Year >= 2010) %>%
    as_tibble() %>%
    select(-OpenInterval) %>% 
    rename_with(tolower) %>% 
    gather(female, male, total, key = sex, value = dts)
  
  chunk_e <- 
    readHMDweb(ct, "Exposures_1x1", hmd_us, hmd_pw) %>%
    # filter(Year >= 2010) %>%
    as_tibble() %>%
    select(-OpenInterval) %>% 
    rename_with(tolower) %>% 
    gather(female, male, total, key = sex, value = exposure)
  
  chunk <- 
    chunk_d %>% 
    left_join(chunk_e) %>% 
    mutate(code = ct)
    
  hmd <- 
    hmd %>%
    bind_rows(chunk)
}

unique(hmd$code)

write_rds(hmd, "data_inter/deaths_exposures_hmd.rds")
