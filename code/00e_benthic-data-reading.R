# This is a script for loading in SBC LTER data from RDS (hopefully saves time)

############################################################-
# 0. source -------------------------------------------------
############################################################-

source(here::here("code", "00a_libraries.R"))

############################################################-
# 1. read RDS files -----------------------------------------
############################################################-

biomass <- readRDS(file = here::here("benthic-data-RDS", "biomass.rds"))
biomass_2022 <- readRDS(file = here::here("benthic-data-RDS", "biomass_2022.rds"))
# percov <- readRDS(file = here::here("benthic-data-RDS", "percov.rds"))
percov_2022 <- readRDS(file = here::here("benthic-data-RDS", "percov_2022.rds"))
swath <- readRDS(file = here::here("benthic-data-RDS", "swath.rds"))
irr <- readRDS(file = here::here("benthic-data-RDS", "irr.rds"))
urchins <- readRDS(file = here::here("benthic-data-RDS", "urchins.rds"))
substrate <- readRDS(file = here::here("benthic-data-RDS", "substrate.rds"))
