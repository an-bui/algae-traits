# This is a script for loading in SBC LTER data from RDS (hopefully saves time)

############################################################-
# 0. source -------------------------------------------------
############################################################-

source(here::here("code", "00a_libraries.R"))

############################################################-
# 1. read RDS files -----------------------------------------
############################################################-

readRDS(biomass, file = here::here("benthic-data-RDS", "biomass.rds"))
readRDS(biomass_2022, file = here::here("benthic-data-RDS", "biomass_2022.rds"))
readRDS(percov, file = here::here("benthic-data-RDS", "percov.rds"))
readRDS(percov_2022, file = here::here("benthic-data-RDS", "percov_2022.rds"))
readRDS(swath, file = here::here("benthic-data-RDS", "swath.rds"))
readRDS(irr, file = here::here("benthic-data-RDS", "irr.rds"))
readRDS(urchins, file = here::here("benthic-data-RDS", "urchins.rds"))
readRDS(substrate, file = here::here("benthic-data-RDS", "substrate.rds"))