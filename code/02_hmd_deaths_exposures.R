# ==============================================================================
# 02. Download death counts and exposures (1x1) from the HMD
#     (used by the APC models of Figures 3 and 4)
#
# Credentials: see header of 01_hmd_life_tables.R
# ==============================================================================
rm(list = ls()); gc()
library(HMDHFDplus)
library(tidyverse)

hmd_us <- Sys.getenv("hmd_us")
hmd_pw <- Sys.getenv("hmd_pw")
if (hmd_us == "" | hmd_pw == "") stop("HMD credentials not found: set hmd_us and hmd_pw (see header of 01_hmd_life_tables.R).")

# Only Spain is modelled in the manuscript; the reference countries are
# included for exploratory work (see code/exploratory/).
# To download every HMD population instead: cds <- getHMDcountries()$CNTRY
cds <- c("ESP", "SWE", "NOR", "NLD", "CHE",   # Spain and the reference countries
         "DNK", "GBRTENW", "FRATNP")         # alternative reference set (Supplementary Figure S1)

dir.create("data_inter", showWarnings = FALSE)

hmd <- tibble()
for(ct in cds){
  cat(paste0(ct, "\n"))
  chunk_d <-
    readHMDweb(ct, "Deaths_1x1", hmd_us, hmd_pw) %>%
    as_tibble() %>%
    select(-OpenInterval) %>%
    rename_with(tolower) %>%
    gather(female, male, total, key = sex, value = dts)

  chunk_e <-
    readHMDweb(ct, "Exposures_1x1", hmd_us, hmd_pw) %>%
    as_tibble() %>%
    select(-OpenInterval) %>%
    rename_with(tolower) %>%
    gather(female, male, total, key = sex, value = exposure)

  chunk <-
    chunk_d %>%
    left_join(chunk_e, by = c("year", "age", "sex")) %>%
    mutate(code = ct)

  hmd <- bind_rows(hmd, chunk)
}

write_rds(hmd, "data_inter/deaths_exposures_hmd.rds")
