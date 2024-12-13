# This is a script for loading in SBC LTER data.

############################################################-
# 0. source -------------------------------------------------
############################################################-

source(here::here("code", "00a_libraries.R"))

############################################################-
# 1. benthic cleaning function ------------------------------
############################################################-

benthic_cleaning_fxn <- function(df) {
  df %>% 
    clean_names() %>% 
    # change to lower case
    mutate_at(c("group", "mobility", "growth_morph", "site"), str_to_lower) %>% 
    # # create a sample_ID for each sampling date at each site
    unite("sample_ID", site, year, transect, remove = FALSE) %>% 
    # only include algae
    filter(group == "algae") %>% 
    # make sure that sp_code for Nienburgia andersoniana isn't NA
    mutate(sp_code = case_when(
      scientific_name == "Nienburgia andersoniana" ~ "Nandersoniana",
      TRUE ~ as.character(as.character(sp_code))
    ))
}

############################################################-
# 2. benthic community data ---------------------------------
############################################################-

# ⊣ a. biomass ----------------------------------------------

biomass <- read_csv(here::here("data", "SBC-LTER-benthics", 
                               "Annual_All_Species_Biomass_at_transect_20240501.csv")) %>% 
  benthic_cleaning_fxn() %>% 
  # replace all -99999 values with NA
  mutate(dry_gm2 = replace(dry_gm2, dry_gm2 < 0, NA),
         wm_gm2 = replace(wm_gm2, wm_gm2 < 0, NA),
         density = replace(density, density < 0, NA)) %>% 
  mutate(date = ymd(date))

saveRDS(biomass, file = here::here("benthic-data-RDS", "biomass.rds"))

# biomass_2022 <- read_csv(here::here("data", "SBC-LTER-benthics", 
#                                     "Annual_All_Species_Biomass_at_transect_20230201.csv")) %>% 
#   benthic_cleaning_fxn() %>% 
#   # replace all -99999 values with NA
#   mutate(dry_gm2 = replace(dry_gm2, dry_gm2 < 0, NA),
#          wm_gm2 = replace(wm_gm2, wm_gm2 < 0, NA),
#          density = replace(density, density < 0, NA),
#          percent_cover = replace(percent_cover, percent_cover < 0, NA)) %>% 
#   mutate(date = ymd(date))
# 
# saveRDS(biomass_2022, file = here::here("benthic-data-RDS", "biomass_2022.rds"))

# ⊣ b. percent cover ----------------------------------------

# percov <- read_csv(here::here("data", "SBC-LTER-benthics", 
#                               "Annual_Cover_All_Years_20211020.csv")) %>% 
#   benthic_cleaning_fxn() %>% 
#   mutate(percent_cover = replace(percent_cover, percent_cover < 0, NA))
# 
# saveRDS(percov, file = here::here("benthic-data-RDS", "percov.rds"))
# 
# percov_2022 <- read_csv(here::here("data", "SBC-LTER-benthics", 
#                                    "Annual_Cover_All_Years_20220809.csv")) %>% 
#   benthic_cleaning_fxn() 
# 
# saveRDS(percov_2022, file = here::here("benthic-data-RDS", "percov_2022.rds"))

# ⊣ c. swath ------------------------------------------------

swath <- read_csv(here::here("data", "SBC-LTER-benthics", 
                             "Annual_Quad_Swath_All_Years_20220809.csv")) %>% 
  benthic_cleaning_fxn() 

saveRDS(swath, file = here::here("benthic-data-RDS", "swath.rds"))


# ⊣ d. LTE --------------------------------------------------

LTE <- read_csv(here::here("data", "SBC-LTER-benthics",
                           "LTE_All_Species_Biomass_at_transect_20220208.csv")) %>% 
  benthic_cleaning_fxn()

saveRDS(LTE, file = here::here("benthic-data-RDS", "LTE.rds"))


############################################################-
# 3. "environmental" data -----------------------------------
############################################################-

# ⊣ a. irradiance -------------------------------------------

# from Castorani et al. 2018: sum irradiance in a day to get total, then average across season to get average daily

irr <- read_csv(here::here("data/SBC-LTER-benthics", "Hourly_Irrandiance_All_Year_20220131.csv")) %>%
  clean_names() %>% 
  # split up dates into year, month, and day
  separate(date_local, into = c("year", "month", "day"), sep = "-", remove = FALSE) %>% 
  # change site and sensor_location to be in lower case
  mutate(across(.cols = c(site, sensor_location), .fns = str_to_lower)) %>% 
  filter(site %in% sites_proposal) %>% 
  filter(!(transect %in% c("MKI", "MKO"))) %>% 
  mutate(across(.cols = c(year, month, day), .fns = as.numeric)) %>% 
  filter(month %in% c(7, 8, 9, 10)) %>% 
  group_by(site, year) %>% 
  summarize(mean_umol = mean(light_umol, na.rm = TRUE),
            se_umol = se(light_umol)) %>% 
  ungroup() %>% 
  unite("sample_ID", site, year, sep = "_")

saveRDS(irr, file = here::here("benthic-data-RDS", "irr.rds"))

# ⊣ b. urchins ----------------------------------------------

# urchin summary
urchins <- read_csv(here::here("data", "SBC-LTER-benthics", 
                               "Annual_All_Species_Biomass_at_transect_20211020.csv")) %>% 
  clean_names() %>% 
  filter(sp_code == "SPL") %>% 
  mutate_at(c("site"), str_to_lower) %>% 
  filter(site %in% sites_proposal) %>% 
  group_by(site, year) %>% 
  # replace all -99999 values with NA
  mutate(dry_gm2 = replace(dry_gm2, dry_gm2 < 0, NA),
         wm_gm2 = replace(wm_gm2, wm_gm2 < 0, NA),
         density = replace(density, density < 0, NA)) %>% 
  summarize(mean_urchins = mean(dry_gm2, na.rm = TRUE),
            se_urchins = se(dry_gm2),
            sum_urchins = sum(dry_gm2)) %>% 
  # create a sample_ID for each sampling date at each site
  unite("sample_ID", site, year, remove = FALSE) %>% 
  ungroup()

saveRDS(urchins, file = here::here("benthic-data-RDS", "urchins.rds"))

# ⊣ c. substrate --------------------------------------------

substrate <- read_csv(here::here("data/SBC-LTER-benthics", "Annual_Substrate_All_Years_20211020.csv")) %>% 
  clean_names() %>% 
  # change to lower case
  mutate_at(c("site"), str_to_lower) %>% 
  filter(site %in% sites_proposal) %>%
  filter(substrate_type == "S") %>% 
  group_by(site, year) %>% 
  summarize(mean_sand = mean(percent_cover, na.rm = TRUE),
            se_sand = se(percent_cover)) %>% 
  # create a sample_ID for each sampling date at each site
  unite("sample_ID", site, year, remove = FALSE) %>% 
  ungroup()

saveRDS(substrate, file = here::here("benthic-data-RDS", "substrate.rds"))



