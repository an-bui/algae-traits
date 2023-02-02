
# 1.  libraries -----------------------------------------------------------

library(tidyverse) # general use
library(googlesheets4) # getting data from google sheets
library(here) # file organization
library(vegan) # ordinations etc
library(corrplot) # correlation plots
library(cluster) # to get gower distance on categorical variables
library(FD) # functional diversity
library(janitor) # cleaning names in data frames
library(lubridate) # dealing with dates
library(gt) # making tables
library(plotly) # for making the species biomass plot
library(calecopal) # colors
library(patchwork) # putting plots together
library(lme4) # GLMERs
library(nlme) # non-linear mixed effect models (gives significance for terms)
library(MuMIn) # lms
library(glmmTMB) # GLMM
library(emmeans) # effect sizes
library(DHARMa) # plotting residuals from `glmmTMB`
library(multcompView) # pairwise comparisons on plots
library(pairwiseAdonis) # pairwise comparisons for permanova
library(ggeffects) # plot model predictions
library(emmeans) # also plot model predictions
library(performance) # checks of collinearity, amongst other things
library(car) # checking VIR
library(cati) # partitioning species composition/intraspecific variation
# library(BiodiversityR) # rank abundance
library(gtsummary)
library(ggnewscale)
library(factoextra)
library(ggridges)


# 2.  getting data from google drive --------------------------------------

# gets the file data from google drive
sheet_id <- "1h2eHoL5kXMRwExt3hjrLavxG0pajHOvA1UfJd24pRk4"

# metadata for subsamples
metadata_sub <- read_sheet(sheet_id, sheet = "00b-metadata_sub") %>% 
  mutate(date_collected = ymd(date_collected)) %>% 
  mutate(year = year(date_collected))

# metadata for individuals 
metadata_ind <- metadata_sub %>% 
  select(specimen_ID, date_collected, year, site, sp_code, lifestage) %>% 
  unique()

# individual maximum height
ind_height <- read_sheet(sheet_id, sheet = "02a-ind_height", na = "NA") 

# thallus length and width
lw <- read_sheet(sheet_id, sheet = "02b-lw", na = "NA") 

# thallus thickness
thickness <- read_sheet(sheet_id, sheet = "03-thickness", na = "NA")

# wet and dry weight
weight <- read_sheet(sheet_id, sheet = "04-weight", na = "NA")

# surface area and perimeter
sa_peri <- read_sheet(sheet_id, sheet = "05a-scans") 

# branching order (still need to clean up)
bra_ord <- read_sheet(sheet_id, sheet = "05b-branching_order", na = "NA") 

# toughness (still need to clean up)
toughness <- read_sheet(sheet_id, sheet = "06-toughness", na = "NA") 

# volume
volume <- read_sheet(sheet_id, sheet = "07-volume", na = "NA") 


# 3.  getting FvFm data ---------------------------------------------------

# creating a vector of file names for .csv
path <- here::here("data/fvfm", c(
  "20210630-MOHK-v2-cleaned.csv",
  "20210709-MOHK-v2-cleaned.csv",
  "20210719-IVEE-v2-cleaned.csv",
  "20210720-CARP-v2-cleaned.csv",
  "20210721-BULL-v2-cleaned.csv",
  "20210722-AQUE-v2-cleaned.csv",
  "20220405-MOHK-cleaned.csv",
  "20220407-AQUE-cleaned.csv",
  "20220421-AQUE-cleaned.csv",
  "20220428-AQUE-cleaned.csv",
  "20220503-AQUE-cleaned.csv",
  "20220719-NAPL-cleaned.csv",
  "20220721-NAPL-cleaned.csv",
  "20220726-BULL-cleaned.csv", 
  "20220811-AQUE-cleaned.csv",
  "20220812-AQUE-cleaned.csv",
  "20220815-MOHK-cleaned.csv",
  "20220818-AQUE-cleaned.csv",
  "20220824-CARP_MOHK-cleaned.csv"))

# read in all .csv files into one data frame
fvfm_raw <- path %>% 
  map_df(~ read_csv(.))

# 4. useful objects -------------------------------------------------------

# vector: most abundant algae
algae_common <- c("PH", "PTCA", # Pterygophora californica 
                  "DL", # Desmarestia ligulata
                  "R", # Rhodymenia californica 
                  "CC", # Chondracanthus corymbiferus 
                  "POLA", # Polyneura latissima 
                  "CYOS", # Stephanocystis osmundacea 
                  "FTHR", # Pterosiphonia dendroidea 
                  "CO", # Corallina officinalis var. chilensis 
                  "LX", # Osmundea spectabilis
                  "GS", # Gracilaria spp. 
                  "GR", # Gelidium robustum
                  "BR", # Halymenia spp.
                  "BO", # Bossiella orbigniana 
                  "FB", # Ectocarpaceae spp. 
                  "BF", # Cryptopleura ruprechtiana 
                  "LAFA", # Laminaria farlowii 
                  "CF", # Callophyllis rhynchocarpa 
                  "DP" # Dictyota spp. 
)

# 11 species from Miller et al. 2012 and from conversation with Bob on 2022-01-18
algae_interest <- c("CYOS", # Stephanocystis osmundacea 
                    "LAFA", # Laminaria farlowii 
                    "MAPY", # Macrocystis pyrifera
                    "PH", "PTCA", # Pterygophora californica 
                    "CF", # Callophyllis flabellulata
                    "CC", # Chondracanthus corymbiferus 
                    "GS", # Gracilaria spp. 
                    "POLA", # Polyneura latissima 
                    "FTHR", # Pterosiphonia dendroidea 
                    "R", # Rhodymenia californica 
                    "EGME", # Egregia menziesii
                    "DL" # Desmarestia ligulata
)

# algae list in proposal
# updated 2023-01-27 with new species
algae_proposal <- c("PH", "PTCA", # Pterygophora californica 
                    "BF", # Cryptopleura ruprechtiana                     
                    "CYOS", # Stephanocystis osmundacea                     
                    "DL", # Desmarestia ligulata                   
                    "CC", # Chondracanthus corymbiferus                     
                    "GS", # Gracilaria spp.                     
                    "CO", # Corallina officinalis var. chilensis 
                    "POLA", # Polyneura latissima 
                    "R", # Rhodymenia californica 
                    "GR", # Gelidium robustum
                    "EH", "EGME", # Egregia menziesii
                    "Nandersoniana", # Nienburgia andersoniana
                    "LH", "LAFA", # Laminaria farlowii
                    "DU", # Dictyopteris undulata
                    "DP", # Dictyota
                    "BO" # Bossiella orbigniana
)

# date
todays_date <- Sys.Date()

# colors
rhodo_col <- "#781416"
ochro_col <- "#CC7540"
chloro_col <- "#6D5A18"

# site name vector
sites <- c("bull", "aque", "ahnd", "napl", "ivee", "golb", 
           "abur", "mohk", "carp", "scdi", "sctw")

sites_proposal <- c("bull", "aque", "napl", "ivee", "mohk", "carp")
sites_proposal_new <- c("aque", "napl", "ivee", "mohk", "carp")

sites_full <- setNames(c("Bullito (BULL)", 
                         "Arroyo Quemado (AQUE)",
                         "Arroyo Hondo (AHND)",
                         "Naples (NAPL)",
                         "Isla Vista (IVEE)",
                         "Goleta Beach (GOLB)",
                         "Arroyo Burro (ABUR)",
                         "Mohawk (MOHK)",
                         "Carpinteria (CARP)",
                         "Diablo Canyon (SCDI)",
                         "Twin Harbors (SCTW)"), sites)

# full names for sites
bull_full <- "Bullito (BULL)"
aque_full <- "Arroyo Quemado (AQUE)"
ahnd_full <- "Arroyo Hondo (AHND)"
napl_full <- "Naples (NAPL)"
ivee_full <- "Isla Vista (IVEE)"
golb_full <- "Goleta Beach (GOLB)"
abur_full <- "Arroyo Burro (ABUR)"
mohk_full <- "Mohawk (MOHK)"
carp_full <- "Carpinteria (CARP)"
scdi_full <- "Diablo Canyon (SCDI)"
sctw_full <- "Twin Harbors (SCTW)"

# gradient palette
gradient_palette <- c("#FFFFFF", "#009BB0")

# function to calculate standard error
se <- function(x,...){
  sd(x, na.rm = TRUE)/sqrt(length(na.omit(x)))
}


# 5.  other files ---------------------------------------------------------

#### * a. traits ####
# coarse trait data
coarse_traits <- read_csv(here::here("data", "algae-traits_literature-search_2022-02-28.csv"))

# intermediate data frame for trait by species matrix
ct_prep <- metadata_sub %>% 
  select(specimen_ID, sp_code, lifestage) %>% 
  left_join(., coarse_traits, by = "sp_code") %>% 
  # select(-subsample_ID) %>% 
  select(specimen_ID, sp_code, scientific_name, 
         growth_form, pigment_type, life_habit, longevity, posture, branching_yn, lifestage) %>%
  unique()

#### * b. benthics ####
benthic_cleaning_fxn <- function(df) {
  df %>% 
    clean_names() %>% 
    # change to lower case
    mutate_at(c("group", "mobility", "growth_morph", "site"), str_to_lower) %>% 
    # # create a sample_ID for each sampling date at each site
    unite("sample_ID", site, year, remove = FALSE) %>% 
    # only include algae
    filter(group == "algae") %>% 
    # make sure that sp_code for Nienburgia andersoniana isn't NA
    mutate(sp_code = case_when(
      scientific_name == "Nienburgia andersoniana" ~ "Nandersoniana",
      TRUE ~ as.character(as.character(sp_code))
    ))
}

#### * c. biomass ####
biomass <- read_csv(here::here("data", "SBC-LTER-benthics", 
                               "Annual_All_Species_Biomass_at_transect_20211020.csv")) %>% 
  benthic_cleaning_fxn() %>% 
  # replace all -99999 values with NA
  mutate(dry_gm2 = replace(dry_gm2, dry_gm2 < 0, NA),
         wm_gm2 = replace(wm_gm2, wm_gm2 < 0, NA),
         density = replace(density, density < 0, NA)) %>% 
  mutate(date = ymd(date))

#### * d. percent cover ####
percov <- read_csv(here::here("data", "SBC-LTER-benthics", 
                              "Annual_Cover_All_Years_20211020.csv")) %>% 
  benthic_cleaning_fxn() %>% 
  mutate(percent_cover = replace(percent_cover, percent_cover < 0, NA))

#### * e. prop_biomass ####
# This is a data frame of species at each site at each transect during each sampling date
# biomass columns: total dry, total wet, percent dry, percent wet
prop_biomass <- biomass %>% 
  # filter(sample_ID == "CARP_2000-09-08") %>% 
  filter(sp_code != "MAPY") %>% 
  filter(sp_code %in% algae_proposal) %>% 
  filter(site %in% sites_proposal) %>% 
  select(sample_ID, site, year, sp_code, dry_gm2, wm_gm2) %>% 
  filter(dry_gm2 > 0) %>% 
  # create a new column where total biomass for the whole site is calculated
  group_by(sample_ID, site, year, sp_code) %>% 
  summarize(total_sp_dry = sum(dry_gm2),
         total_sp_wet = sum(wm_gm2)) %>%
  ungroup() %>% 
  # create a new column where total biomass for the species for the whole site is calculated
  group_by(site, year) %>%
  mutate(total_dry = sum(total_sp_dry),
         total_wet = sum(total_sp_wet)) %>%
  ungroup() %>%
  # create a new column where percent of total biomass per survey for each species is calculated
  mutate(percent_sp_dry = total_sp_dry/total_dry,
         percent_sp_wet = total_sp_wet/total_wet) %>% 
  # replace NaNs with 0 (which happens when there's nothing in the survey)
  # mutate_at(vars("percent_sp_dry", "percent_sp_wet"), list(~ ifelse(. = NaN, 0, .))) %>%
  # double check to make sure the percents worked out...
  select(sample_ID, site, year, sp_code, 
         # values that need to be unique: 
         total_dry, total_wet, percent_sp_dry, percent_sp_wet) %>% 
  # unique() %>% 
  group_by(site, year) %>% 
  mutate(total_percent_dry = sum(percent_sp_dry),
         total_percent_wet = sum(percent_sp_wet)) %>% 
  # join with coarse traits data frame, which includes taxonomy
  left_join(., coarse_traits, by = "sp_code") %>% 
  ungroup()

av_biomass <- biomass %>% 
  # filter(sample_ID == "CARP_2000-09-08") %>% 
  filter(sp_code %in% algae_proposal) %>% 
  filter(site %in% sites_proposal) %>% 
  group_by(site, year, sp_code) %>% 
  summarize(mean_dry = mean(dry_gm2, na.rm = TRUE)) %>% 
  ungroup() %>% 
  unite("sample_ID", site, year, remove = FALSE) %>% 
  # join with coarse traits data frame, which includes taxonomy
  left_join(., coarse_traits, by = "sp_code")

biomass_transect <- biomass %>% 
  dplyr::select(-sample_ID) %>% 
  unite("sample_ID", site, year, remove = FALSE) %>% 
  unite("transect_ID", site, year, transect, remove = FALSE) %>% 
  filter(sp_code %in% algae_proposal) %>% 
  filter(site %in% sites_proposal) %>% 
  # join with coarse traits data frame, which includes taxonomy
  left_join(., coarse_traits, by = "sp_code") 

#### * f. irradiance ####

# * f. irradiance
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

#### * g. community matrix and metadata ####

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

urchin_summary <- urchins %>%
  select(sample_ID, mean_urchins, se_urchins, sum_urchins)

urchin_transect_summary <- read_csv(here::here("data", "SBC-LTER-benthics", 
  "Annual_All_Species_Biomass_at_transect_20211020.csv")) %>% 
  clean_names() %>% 
  filter(sp_code == "SPL") %>% 
  mutate_at(c("site"), str_to_lower) %>% 
  filter(site %in% sites_proposal) %>% 
  group_by(site, year, transect) %>% 
  summarize(mean_urchins = mean(dry_gm2, na.rm = TRUE),
            se_urchins = se(dry_gm2),
            sum_urchins = sum(dry_gm2)) %>% 
  # create a sample_ID for each sampling date at each site
  # unite("sample_ID", site, year, remove = FALSE) %>% 
  unite("transect_ID", site, year, transect, remove = TRUE) %>% 
  ungroup()

# sand summary
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

substrate_transect_summary <- read_csv(here::here("data/SBC-LTER-benthics", "Annual_Substrate_All_Years_20211020.csv")) %>% 
  clean_names() %>% 
  # change to lower case
  mutate_at(c("site"), str_to_lower) %>% 
  filter(site %in% sites_proposal) %>%
  filter(!(common_name %in% c("Sand", "Shell debris", "Shallow Sand"))) %>% 
  group_by(site, year, transect) %>% 
  summarize(mean_subs = mean(percent_cover, na.rm = TRUE),
            se_subs = se(percent_cover),
            max_subs = max(percent_cover)) %>% 
  # create a sample_ID for each sampling date at each site
  unite("transect_ID", site, year, transect, remove = TRUE) %>% 
  ungroup()

# kelp summary
kelp_biomass_summary <- biomass %>% 
  filter(sp_code == "MAPY") %>% 
  filter(site %in% sites_proposal) %>% 
  # replace all -99999 values with 0
  mutate(wm_gm2 = replace(wm_gm2, wm_gm2 < 0, NA),
         dry_gm2 = replace(dry_gm2, dry_gm2 < 0, NA)) %>% 
  dplyr::select(sample_ID, dry_gm2, wm_gm2) %>% 
  rename(kelp_dry = dry_gm2,
         kelp_wet = wm_gm2) %>% 
  group_by(sample_ID) %>%
  summarize(mean_kelp = mean(kelp_dry, na.rm = TRUE),
            se_kelp = se(kelp_dry),
            sum_kelp = sum(kelp_dry)) %>%
  ungroup()

kelp_biomass_transect_summary <- biomass %>% 
  dplyr::select(-sample_ID) %>% 
  unite("sample_ID", site, year, remove = FALSE) %>% 
  unite("transect_ID", site, year, transect, remove = FALSE) %>% 
  filter(sp_code == "MAPY") %>% 
  filter(site %in% sites_proposal) %>% 
  # replace all -99999 values with 0
  mutate(wm_gm2 = replace(wm_gm2, wm_gm2 < 0, NA),
         dry_gm2 = replace(dry_gm2, dry_gm2 < 0, NA)) %>% 
  rename(kelp_dry = dry_gm2,
         kelp_wet = wm_gm2) %>% 
  group_by(site, year, transect) %>%
  summarize(mean_kelp = mean(kelp_dry, na.rm = TRUE),
            se_kelp = se(kelp_dry),
            sum_kelp = sum(kelp_dry)) %>% 
  ungroup() %>% 
  unite("transect_ID", site, year, transect, remove = TRUE)

# site by species matrix
community_matrix <- prop_biomass %>% 
  # filter out sites of interest
  # filter(site %in% c("bull", "aque", "napl", "ivee", "mohk", "carp")) %>% 
  # select percent biomass
  dplyr::select(sample_ID, sp_code, percent_sp_dry) %>% 
  pivot_wider(names_from = sp_code, values_from = percent_sp_dry) %>% 
  # take out surveys with 0 observations
  # rowwise() %>% 
  # mutate(sum = sum(c_across(BF:R))) %>% 
  # ungroup() %>% 
  # select(-sum) %>% 
  column_to_rownames("sample_ID") %>% 
  mutate_all(~replace(., is.na(.), 0))

community_matrix_av <- av_biomass %>% 
  dplyr::select(sample_ID, sp_code, mean_dry) %>% 
  pivot_wider(names_from = sp_code, values_from = mean_dry) %>% 
  column_to_rownames("sample_ID") %>% 
  mutate_all(~replace(., is.na(.), 0))

community_matrix_transect <- biomass_transect %>% 
  dplyr::select(transect_ID, sp_code, dry_gm2) %>% 
  pivot_wider(names_from = sp_code, values_from = dry_gm2) %>% 
  column_to_rownames("transect_ID") %>% 
  mutate_all(~replace(., is.na(.), 0))

# metadata for plotting
community_metadata <- biomass %>% 
  # filter out sites of interest
  filter(site %in% sites_proposal) %>% 
  dplyr::select(sample_ID, site, year) %>% 
  unique() %>% 
  filter(sample_ID %in% rownames(community_matrix)) %>% 
  left_join(., urchin_summary, by = "sample_ID") %>% 
  left_join(., substrate, by = "sample_ID") %>% 
  left_join(., kelp_biomass_summary, by = "sample_ID") %>% 
  left_join(., irr, by = c("sample_ID"))

community_metadata_transect <- biomass %>% 
  # filter out sites of interest
  filter(site %in% sites_proposal) %>% 
  dplyr::select(site, year, transect) %>% 
  unique() %>% 
  unite("sample_ID", site, year, remove = FALSE) %>% 
  unite("transect_ID", site, year, transect, remove = FALSE) %>% 
  filter(transect_ID %in% rownames(community_matrix_transect)) %>% 
  left_join(., urchin_transect_summary, by = "transect_ID") %>% 
  left_join(., substrate_transect_summary, by = "transect_ID") %>% 
  left_join(., kelp_biomass_transect_summary, by = "transect_ID") %>% 
  left_join(., irr, by = c("sample_ID"))

community_metadata_transect_sub <- community_metadata_transect %>% 
  drop_na(mean_umol)

community_matrix_transect_sub <- biomass_transect %>% 
  dplyr::select(transect_ID, sp_code, dry_gm2) %>% 
  filter(transect_ID %in% community_metadata_transect_sub$transect_ID) %>% 
  pivot_wider(names_from = sp_code, values_from = dry_gm2) %>% 
  # take out surveys with 0 observations
  rowwise() %>%
  mutate(sum = sum(c_across(BF:R))) %>%
  filter(sum > 0) %>% 
  dplyr::select(-sum) %>% 
  ungroup() %>% 
  column_to_rownames("transect_ID") %>% 
  mutate_all(~replace(., is.na(.), 0))

community_metadata_transect_sub <- community_metadata_transect %>% 
  filter(transect_ID %in% rownames(community_matrix_transect_sub))


# check to make sure rows match up
# rownames(community_matrix) == community_metadata$sample_ID

# community matrix in presence-absence form
community_presabs <- community_matrix_transect %>% 
  mutate(across(.cols = everything(), ~replace(., . > 0, 1)))



# 6. trait cleaning/summarizing -------------------------------------------

# This section puts traits into a format that can be used in `cati` and ordinate.

#### * a. Fv/Fm ####

# fvfm for each subsample
fvfm_sub <- fvfm_raw %>% 
  rename('fvfm_meas' = '1:Fv/Fm') %>% 
  # fill in '1:Fv/Fm' with 1:Y(II) when there is none
  # only an issue with 20220405-MOHK samples
  mutate(fvfm_meas = case_when(
    !is.na((fvfm_meas)) ~ `1:Y (II)`,
    TRUE ~ fvfm_meas
  )) %>% 
  # take out observations with "-" in the measurement
  filter(fvfm_meas != "-") %>% 
  # only use FvFm, drop the NAs, change the column name
  select(specimen_ID, subsample_ID, fvfm_meas) %>% 
  drop_na(specimen_ID) %>% 
  # make sure that fvfm_meas is numeric
  mutate(fvfm_meas = as.numeric(fvfm_meas)) %>% 
  # calculate mean, variance, standard deviation
  group_by(specimen_ID, subsample_ID) %>% 
  summarize(fvfm_mean = mean(fvfm_meas, na.rm = TRUE),
         fvfm_se = se(fvfm_meas)) %>% 
  ungroup() %>% 
  # drop NAs in subsample_ID: 20210719-IVEE-026 and 20210630-MOHK-007
  drop_na(subsample_ID) %>% 
  select(-specimen_ID)

# fvfm for each individual
fvfm_ind <- fvfm_raw %>% 
  rename('fvfm_meas' = '1:Fv/Fm') %>% 
  # fill in '1:Fv/Fm' with 1:Y(III) when there is none
  # only an issue with 20220405-MOHK samples
  mutate(fvfm_meas = case_when(
    !is.na((fvfm_meas)) ~ `1:Y (II)`,
    TRUE ~ fvfm_meas
  )) %>% 
  # take out observations with "-" in the measurement
  filter(fvfm_meas != "-")%>% 
  # only use FvFm, drop the NAs, change the column name
  select(specimen_ID, subsample_ID, 'fvfm_meas') %>% 
  drop_na(specimen_ID) %>%
  # make sure that fvfm_meas is numeric
  mutate(fvfm_meas = as.numeric(fvfm_meas)) %>% 
  # calculate mean, variance, standard deviation
  group_by(specimen_ID) %>% 
  summarize(fvfm_mean = mean(fvfm_meas, na.rm = TRUE),
            fvfm_se = se(fvfm_meas)) %>% 
  ungroup() 

# as a note: samples from 20210623 don't have FvFm, and one sample from IVEE doesn't either

#### * b. thickness ####

# thickness per each subsample
thickness_sub <- thickness %>% 
  pivot_longer(cols = thickness_01:thickness_10, names_to = "measurement_n", values_to = "thickness_mm") %>% 
  group_by(subsample_ID) %>% 
  summarize(thickness_mm_mean = mean(thickness_mm, na.rm = TRUE),
            thickness_mm_se = se(thickness_mm))

# thickness per individual
thickness_ind <- thickness %>% 
  pivot_longer(cols = thickness_01:thickness_10, names_to = "measurement_n", values_to = "thickness_mm") %>% 
  group_by(specimen_ID) %>% 
  summarize(thickness_mm_mean = mean(thickness_mm, na.rm = TRUE),
            thickness_mm_se = se(thickness_mm))

#### * c. weights ####

# weight of each subsample
weight_sub <- weight %>%
  select(specimen_ID, subsample_ID, weight_wet_mg, weight_wet_g, weight_dry_mg, weight_dry_g) %>%
  # calculate dry matter content
  mutate(dmc = weight_dry_mg/weight_wet_mg) %>%
  # take out specimen_ID column
  select(-specimen_ID)

# weight of each individual
weight_ind <- weight %>% 
  group_by(specimen_ID) %>% 
  summarize(total_wet = sum(weight_wet_mg, na.rm = TRUE),
            total_dry = sum(weight_dry_mg, na.rm = TRUE)) %>% 
  ungroup()

#### * d. volume ####

# volume of each subsample
volume_sub <- volume %>%
  select(subsample_ID, volume_total_mL)

# volume of each individual
volume_ind <- volume %>% 
  group_by(specimen_ID) %>% 
  summarize(total_volume = sum(volume_total_mL, na.rm = TRUE))

#### * e. SA, perimeter ####

# surface area and perimeter of each subsample
sa_peri_sub <- sa_peri %>% 
  select(specimen_ID, subsample_ID, area_total, peri_total) %>% 
  group_by(subsample_ID) %>% 
  summarize(area_total = sum(area_total),
            peri_total = sum(peri_total)) 

#### * f. max height and width ####

# thallus length and width
lw_sub <- lw %>% 
  select(-specimen_ID, -5)

#### * g. branching order ####

# bra_ord_summary <- bra_ord %>% 
#   pivot_longer(cols = bo_01:bo_05, names_to = "measurement_n", values_to = "bra_ord_count") %>% 
#   # calculate mean branching order for each sample
#   group_by(specimen_ID) %>% 
#   summarize(mean = mean(bra_ord_count, na.rm = TRUE)) %>% 
#   # join with metadata_subsamples
#   left_join(., metadata_subsamples, by = "specimen_ID") %>% 
#   # join with coarse_traits
#   left_join(., coarse_traits, by = "sp_code") %>% 
#   drop_na(sp_code)

#### * h. toughness ####

# toughness_summary <- toughness %>% 
#   pivot_longer(cols = 2:11, names_to = "measurement_n", values_to = "toughness_kgcm2") %>% 
#   # calculate mean toughness for each sample
#   group_by(specimen_ID) %>% 
#   summarize(mean = mean(toughness_kgcm2, na.rm = TRUE)) %>% 
#   # join with metadata_subsamples
#   left_join(., metadata_subsamples, by = "specimen_ID") %>% 
#   # join with coarse_traits
#   left_join(., coarse_traits, by = "sp_code") %>% 
#   drop_na(sp_code)

#### * i. trait by sample matrix ####

# leaf traits only
leaf_traits <- metadata_sub %>% 
  filter(type %in% c("whole", "thallus")) %>% 
  left_join(., fvfm_sub, by = "subsample_ID") %>% 
  left_join(., thickness_sub, by = "subsample_ID") %>% 
  left_join(., weight_sub, by = "subsample_ID") %>% 
  left_join(., volume_sub, by = "subsample_ID") %>% 
  left_join(., lw_sub, by = "subsample_ID") %>% 
  left_join(., sa_peri_sub, by = "subsample_ID") %>% 
  # ratios
  mutate(sap_ratio = area_total/peri_total,
         sav_ratio = area_total/volume_total_mL,
         sta_mm_mg = area_total/weight_dry_mg)

# average leaf values
av_leaf_values <- leaf_traits %>% 
  group_by(specimen_ID) %>% 
  summarize(sta_mean = mean(sta_mm_mg, na.rm = TRUE),
            sta_se = se(sta_mm_mg),
            tdmc_mean = mean(dmc, na.rm = TRUE),
            tdmc_se = se(dmc),
            sav_mean = mean(sav_ratio, na.rm = TRUE),
            sav_se = se(sav_ratio),
            sap_mean = mean(sap_ratio, na.rm = TRUE),
            sap_se = se(sap_ratio),
            frond_length_mean = mean(length, na.rm = TRUE),
            frond_length_se = se(length),
            frond_width_mean = mean(width, na.rm = TRUE), 
            frond_width_se = se(width))

# all traits for each individual with mean taken for each trait
ind_traits <- ct_prep %>% 
  left_join(., ind_height, by = "specimen_ID") %>% 
  left_join(., thickness_ind, by = "specimen_ID") %>% 
  left_join(., av_leaf_values, by = "specimen_ID") %>% 
  left_join(., fvfm_ind, by = "specimen_ID") %>% 
  left_join(., weight_ind, by = "specimen_ID") %>% 
  left_join(., volume_ind, by = "specimen_ID") %>% 
  left_join(., (metadata_ind %>% select(specimen_ID, date_collected, site)), by = "specimen_ID") %>% 
  filter(!(sp_code == "PTCA" & lifestage == "recruit")) %>% 
  mutate(date_collected = ymd(date_collected)) %>% 
  mutate(year = year(date_collected))

#### * j. trait by species matrix ####

# trait values averaged across sites
tbspp_matrix <- metadata_ind %>% 
  select(specimen_ID) %>% 
  unique() %>% 
  left_join(., ind_height, by = "specimen_ID") %>% 
  left_join(., thickness_ind, by = "specimen_ID") %>% 
  left_join(., av_leaf_values, by = "specimen_ID") %>% 
  left_join(., fvfm_ind, by = "specimen_ID") %>% 
  left_join(., weight_ind, by = "specimen_ID") %>% 
  left_join(., volume_ind, by = "specimen_ID") %>% 
  left_join(., metadata_ind, by = "specimen_ID") %>% 
  filter(!(sp_code == "PTCA" & lifestage == "recruit")) %>% 
  select(specimen_ID, sp_code,
         fvfm_mean, thickness_mm_mean, total_dry, tdmc_mean,
         total_volume, sta_mean, sap_mean, sav_mean, maximum_height) %>% 
  # calculate mean trait value for each species
  group_by(sp_code) %>% 
  summarize_at(vars(fvfm_mean:maximum_height), mean, na.rm = TRUE) %>% 
  ungroup() %>% 
  drop_na(sp_code) %>% 
  # join with coarse traits
  left_join(., coarse_traits, by = "sp_code") %>% 
  select(sp_code:taxon_phylum, growth_form, pigment_type, life_habit, longevity, posture, branching_yn) %>% 
  filter(sp_code != "MAPY") %>% 
  filter(sp_code %in% algae_proposal) %>% 
  column_to_rownames("sp_code") %>% 
  # replace NaNs with NAs
  mutate_all(na_if, "NaN") 
  

# tbspp_matrix <- full_join(fvfm_prep, thickness_prep, by = "specimen_ID") %>% 
#   full_join(., weight_prep, by = "specimen_ID") %>% 
#   full_join(., volume_prep, by = "specimen_ID") %>% 
#   full_join(., sa_peri_prep, by = "specimen_ID") %>% 
#   full_join(., hw_prep, by = "specimen_ID") %>% 
#   full_join(., sta_prep, by = "specimen_ID") %>% 
#   full_join(., ct_prep, by = "specimen_ID") %>% 
#   select(specimen_ID, sp_code,
#          fvfm, thickness_mm, weight_dry_g, tdmc, volume_mL,
#          area_total, peri_total, sap_ratio, max_height_cm, max_width_cm, sta_mm_mg) %>% 
#   # calculate mean trait value for each species
#   group_by(sp_code) %>% 
#   summarize_at(vars(fvfm:sta_mm_mg), mean, na.rm = TRUE) %>% 
#   ungroup() %>% 
#   drop_na(sp_code) %>% 
#   # join with coarse_traits 
#   left_join(., coarse_traits, by = "sp_code") %>% 
#   select(sp_code, 
#          fvfm, thickness_mm.x, weight_dry_g, tdmc, volume_mL,
#          area_total, peri_total, sap_ratio, max_height_cm, max_width_cm, sta_mm_mg, 
#          # categorical
#          growth_form, pigment_type, life_habit, longevity, posture, branching_yn) %>% 
#   rename(thickness_mm = thickness_mm.x) %>% 
#   filter(sp_code != "MAPY") %>% 
#   column_to_rownames("sp_code") %>% 
#   # replace NaNs with NAs
#   mutate_all(na_if, "NaN")

#### * k. trait by species matrix function ####

specif_tbspp <- function(site_choice) {
  df <- metadata_ind %>% 
    select(specimen_ID) %>% 
    unique() %>% 
    left_join(., ind_height, by = "specimen_ID") %>% 
    left_join(., thickness_ind, by = "specimen_ID") %>% 
    left_join(., av_leaf_values, by = "specimen_ID") %>% 
    left_join(., fvfm_ind, by = "specimen_ID") %>% 
    left_join(., weight_ind, by = "specimen_ID") %>% 
    left_join(., volume_ind, by = "specimen_ID") %>% 
    left_join(., metadata_ind, by = "specimen_ID") %>% 
    filter(!(sp_code == "PTCA" & lifestage == "recruit")) %>% 
    filter(site == {{ site_choice }}) %>% 
    select(specimen_ID, sp_code,
           fvfm_mean, thickness_mm_mean, total_dry, tdmc_mean,
           total_volume, sta_mean, sap_mean, sav_mean, maximum_height) %>% 
    filter(sp_code != "MAPY") %>% 
    filter(sp_code %in% algae_proposal) %>% 
    # calculate mean trait value for each species
    group_by(sp_code) %>% 
    summarize_at(vars(fvfm_mean:maximum_height), mean, na.rm = TRUE) %>% 
    ungroup() %>% 
    drop_na(sp_code) %>% 
    column_to_rownames("sp_code") %>% 
    # replace NaNs with NAs
    mutate_all(na_if, "NaN")
}

# specif_tbspp <- function(site_choice) {
#   df <- full_join(fvfm_prep, thickness_prep, by = "specimen_ID") %>% 
#     full_join(., weight_prep, by = "specimen_ID") %>% 
#     full_join(., volume_prep, by = "specimen_ID") %>% 
#     full_join(., sa_peri_prep, by = "specimen_ID") %>% 
#     full_join(., hw_prep, by = "specimen_ID") %>% 
#     full_join(., sta_prep, by = "specimen_ID") %>% 
#     full_join(., ct_prep, by = "specimen_ID") %>% 
#     filter(site == {{ site_choice }}) %>% 
#   select(specimen_ID, sp_code,
#          fvfm, thickness_mm, weight_dry_g, tdmc, volume_mL,
#          area_total, peri_total, sap_ratio, max_height_cm, max_width_cm, sta_mm_mg) %>% 
#     # calculate mean trait value for each species
#     group_by(sp_code) %>% 
#     summarize_at(vars(fvfm:sta_mm_mg), mean, na.rm = TRUE) %>% 
#     ungroup() %>% 
#     drop_na(sp_code) %>% 
#     # join with coarse_traits 
#     left_join(., coarse_traits, by = "sp_code") %>% 
#     select(sp_code, 
#            fvfm, thickness_mm.x, weight_dry_g, tdmc, volume_mL,
#            area_total, peri_total, sap_ratio, max_height_cm, max_width_cm, sta_mm_mg, 
#            # categorical
#            growth_form, pigment_type, life_habit, longevity, posture, branching_yn) %>% 
#     rename(thickness_mm = thickness_mm.x) %>% 
#     filter(sp_code != "MAPY") %>% 
#     column_to_rownames("sp_code") %>% 
#     # replace NaNs with NAs
#     mutate_all(na_if, "NaN")
#   
#   return(df)
# }

# note: double check this worked by looking at one of the summary data frames
# and grouping by site and species code, then calculating mean (it does)

#### * l. community matrix function ####

specif_commat <- function(site_choice, specific_tbspp) {
  
  mat <- biomass_transect %>% 
    # filter out sites of interest
    filter(site == {{ site_choice }}) %>% 
    dplyr::select(transect_ID, sp_code, dry_gm2) %>% 
    pivot_wider(names_from = sp_code, values_from = dry_gm2) %>% 
    # column_to_rownames("transect_ID") %>% 
    mutate_all(~replace(., is.na(.), 0)) %>% 
    # take out surveys with 0 observations
    rowwise() %>% 
    mutate(sum = sum(c_across(2:13))) %>% 
    ungroup() %>% 
    filter(sum > 0) %>% 
    select(-sum) %>% 
    column_to_rownames("transect_ID") %>% 
    select(c(rownames(specific_tbspp))) %>% 
    # take out 0 sites
    mutate(total = rowSums(across(where(is.numeric)))) %>% 
    filter(total > 0) %>% 
    select(-total) %>% 
    as.matrix()  
  
  return(mat)
  
}


