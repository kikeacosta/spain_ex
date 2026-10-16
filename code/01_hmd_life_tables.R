# ==============================================================================
# 01. Download period life expectancy at birth and complete (1x1) life tables
#     from the Human Mortality Database (HMD, www.mortality.org)
#
# HMD data are not redistributed with this repository. You need a (free) HMD
# account; store your credentials as environment variables, e.g. in ~/.Renviron:
#   hmd_us = "your_email"
#   hmd_pw = "your_password"
# ==============================================================================
rm(list = ls()); gc()
library(HMDHFDplus)
library(tidyverse)

hmd_us <- Sys.getenv("hmd_us")
hmd_pw <- Sys.getenv("hmd_pw")
if (hmd_us == "" | hmd_pw == "") stop("HMD credentials not found: set hmd_us and hmd_pw (see header).")

# Spain and the four reference countries used in the manuscript.
# To download every HMD population instead: cds <- getHMDcountries()$CNTRY
cds <- c("ESP", "SWE", "NOR", "NLD", "CHE",   # Spain and the reference countries
         "DNK", "GBRTENW", "FRATNP")         # alternative reference set (Supplementary Figure S1)

dir.create("data_inter", showWarnings = FALSE)

# life expectancy at birth by period (Figure 1) ---------------------------------
hmd_e0 <- tibble()
for(ct in cds){
  cat(paste0(ct, "\n"))
  chunk <-
    readHMDweb(ct, "E0per", hmd_us, hmd_pw) %>%
    as_tibble() %>%
    mutate(Code = ct)
  hmd_e0 <- bind_rows(hmd_e0, chunk)
}
write_rds(hmd_e0, "data_inter/e0_hmd.rds")

# complete life tables by sex (Figure 5) ----------------------------------------
read_lt <- function(ct, item, sx){
  readHMDweb(ct, item, hmd_us, hmd_pw) %>%
    filter(Year >= 1900) %>%
    as_tibble() %>%
    mutate(code = ct, sex = sx) %>%
    rename(year = Year, age = Age)
}

hmd_lt <- tibble()
for(ct in cds){
  cat(paste0(ct, "\n"))
  hmd_lt <-
    bind_rows(hmd_lt,
              read_lt(ct, "bltper_1x1", "t"),
              read_lt(ct, "fltper_1x1", "f"),
              read_lt(ct, "mltper_1x1", "m"))
}
write_rds(hmd_lt, "data_inter/ltables_hmd.rds")

