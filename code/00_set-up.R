
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
library(MuMIn) # lms
library(glmmTMB) # GLMM
library(emmeans) # effect sizes
library(DHARMa) # plotting residuals from `glmmTMB`
library(multcompView) # pairwise comparisons on plots
library(pairwiseAdonis) # pairwise comparisons for permanova
library(ggeffects) # plot model predictions

# 2.  getting data from google drive --------------------------------------

# gets the file data from google drive
sheet_id <- "1h2eHoL5kXMRwExt3hjrLavxG0pajHOvA1UfJd24pRk4"

metadata <- read_sheet(sheet_id, sheet = "00b-metadata") %>% 
  filter(site != "Campus Point")

metadata_subsamples <- read_sheet(sheet_id, sheet = "00c-metadata-subsamples") %>% 
  filter(site != "Campus Point")

hw <- read_sheet(sheet_id, sheet = "02-hw", na = "NA") 

thickness <- read_sheet(sheet_id, sheet = "03-thickness", na = "NA") 

weight <- read_sheet(sheet_id, sheet = "04-weight", na = "NA") 

sa_peri <- read_sheet(sheet_id, sheet = "05a-scans") 

bra_ord <- read_sheet(sheet_id, sheet = "05b-branching_order", na = "NA") 

toughness <- read_sheet(sheet_id, sheet = "06-toughness", na = "NA") 

volume <- read_sheet(sheet_id, sheet = "07-volume", na = "NA") 


# 3.  getting FvFm data ---------------------------------------------------

# creating a vector of file names for .csv
path <- here::here("data", c(
  "20210630-MOHK-v2-cleaned.csv",
  "20210709-MOHK-v2-cleaned.csv",
  "20210719-IVEE-v2-cleaned.csv",
  "20210720-CARP-v2-cleaned.csv",
  "20210721-BULL-v2-cleaned.csv",
  "20210722-AQUE-v2-cleaned.csv"))

# read in all .csv files into one data frame
fvfm_raw <- path %>% 
  map_df(~ read_csv(.))


# 4.  other files ---------------------------------------------------------

#### * a. traits ####
# coarse trait data
coarse_traits <- read_csv(here::here("data", "algae-traits_literature-search_2022-02-28.csv"))

## b. benthics
benthic_cleaning_fxn <- function(df) {
  df %>% 
    clean_names() %>% 
    # create a sample_ID for each sampling date at each site
    unite("sample_ID", site, date, transect, remove = FALSE) %>% 
    # change to lower case
    mutate_at(c("group", "mobility", "growth_morph", "site"), str_to_lower) %>% 
    # only include algae
    filter(group == "algae") %>% 
    # make sure that sp_code for Nienburgia andersoniana isn't NA
    mutate(sp_code = case_when(
      scientific_name == "Nienburgia andersoniana" ~ "Nandersoniana",
      TRUE ~ as.character(as.character(sp_code))
    ))
}

#### * b. biomass ####
biomass <- read_csv(here::here("data", "SBC-LTER-benthics", 
                               "Annual_All_Species_Biomass_at_transect_20210108.csv")) %>% 
  benthic_cleaning_fxn() %>% 
  # replace all -99999 values with 0
  mutate(dry_gm2 = replace(dry_gm2, dry_gm2 < 0, 0),
         wm_gm2 = replace(wm_gm2, wm_gm2 < 0, 0),
         density = replace(density, density < 0, 0)) %>% 
  mutate(date = ymd(date))

#### * c. percent cover ####
percov <- read_csv(here::here("data", "SBC-LTER-benthics", 
                              "Annual_Cover_All_Years_20210108.csv")) %>% 
  benthic_cleaning_fxn()

#### * d. prop_biomass ####
# This is a data frame of species at each site at each transect during each sampling date
# biomass columns: total dry, total wet, percent dry, percent wet
prop_biomass <- biomass %>% 
  filter(sp_code != "MAPY") %>% 
  select(sample_ID, site, year, date, transect, sp_code, dry_gm2, wm_gm2) %>% 
  # create a new column where total biomass for the whole survey is calculated
  group_by(sample_ID, year) %>% 
  mutate(total_dry = sum(dry_gm2),
         total_wet = sum(wm_gm2)) %>% 
  ungroup() %>% 
  # create a new column where total biomass for the species for the whole survey is calculated
  # group_by(sample_ID, sp_code) %>% 
  # mutate(total_sp_dry = sum(dry_gm2),
  #        total_sp_wet = sum(wm_gm2)) %>% 
  # ungroup() %>% 
  # create a new column where percent of total biomass per survey for each species is calculated
  mutate(percent_sp_dry = dry_gm2/total_dry,
         percent_sp_wet = wm_gm2/total_wet) %>% 
  # replace NaNs with 0 (which happens when there's nothing in the survey)
  # mutate_at(vars("percent_sp_dry", "percent_sp_wet"), list(~ ifelse(. = NaN, 0, .))) %>% 
  # double check to make sure the percents worked out...
  select(sample_ID, site, year, date, transect, sp_code, 
         # values that need to be unique: 
         total_dry, total_wet, dry_gm2, wm_gm2, percent_sp_dry, percent_sp_wet) %>% 
  unique() %>% 
  group_by(sample_ID) %>% 
  mutate(total_percent_dry = sum(percent_sp_dry),
         total_percent_wet = sum(percent_sp_wet)) %>% 
  # join with coarse traits data frame, which includes taxonomy
  left_join(., coarse_traits, by = "sp_code")

#### * e. community matrix and metadata ####

# urchin summary
urchins <- read_csv(here::here("data", "SBC-LTER-benthics", 
                               "Annual_All_Species_Biomass_at_transect_20210108.csv")) %>% 
  clean_names() %>% 
  # create a sample_ID for each sampling date at each site
  unite("sample_ID", site, date, transect, remove = FALSE) %>% 
  filter(sp_code == "SPL") %>% 
  mutate_at(c("site"), str_to_lower) %>% 
  filter(!(site %in% c("scdi", "sctw")))

urchin_summary <- urchins %>%
  select(sample_ID, density) %>% 
  rename(urchin_density = density)

# sand summary
substrate <- read_csv(here::here("data/SBC-LTER-benthics", "Annual_Substrate_All_Years_20211020.csv")) %>% 
  clean_names() %>% 
  # create a sample_ID for each sampling date at each site
  unite("sample_ID", site, date, transect, remove = FALSE) %>% 
  # change to lower case
  mutate_at(c("site"), str_to_lower) %>% 
  filter(!(site %in% c("sctw", "scdi")))

substrate_summary <- substrate %>%
  filter(common_name == "Sand") %>%
  select(sample_ID, percent_cover) %>% 
  rename(sand_pc = percent_cover) %>% 
  group_by(sample_ID) %>% 
  summarize(mean_sand = mean(sand_pc, na.omit = TRUE))

# kelp summary
kelp_biomass_summary <- biomass %>% 
  filter(sp_code == "MAPY") %>% 
  # replace all -99999 values with 0
  mutate(wm_gm2 = replace(wm_gm2, wm_gm2 < 0, 0)) %>% 
  select(sample_ID, dry_gm2, wm_gm2) %>% 
  rename(kelp_dry = dry_gm2,
         kelp_wet = wm_gm2)
  # group_by(sample_ID) %>% 
  # summarize(mean_kelp = mean(dry_gm2, na.rm = TRUE)) %>% 
  # ungroup()

# site by species matrix
community_matrix <- prop_biomass %>% 
  # filter out sites of interest
  filter(site %in% c("aque", "napl", "ivee", "mohk", "carp")) %>% 
  select(sample_ID, sp_code, dry_gm2) %>% 
  pivot_wider(names_from = "sp_code", values_from = "dry_gm2") %>% 
  # take out surveys with 0 observations
  rowwise() %>% 
  mutate(sum = sum(c_across(1:50))) %>% 
  ungroup() %>% 
  filter(sum > 0) %>% # only 2: CARP_2000-09-08 and NAPL_2000-09-18
  filter(sample_ID != "AQUE_2019-07-23") %>% 
  select(-sum) %>% 
  column_to_rownames("sample_ID")

# metadata for plotting
community_metadata <- biomass %>% 
  # filter out sites of interest
  filter(site %in% c("aque", "napl", "ivee", "mohk", "carp")) %>% 
  # filter out SCI sites and AHND
  # filter(!(site %in% c("sctw", "scdi", "ahnd"))) %>% 
  select(sample_ID, site, year, transect) %>% 
  unique() %>% 
  filter(sample_ID %in% rownames(community_matrix)) %>% 
  left_join(., urchin_summary, by = "sample_ID") %>% 
  left_join(., substrate_summary, by = "sample_ID") %>% 
  left_join(., kelp_biomass_summary, by = "sample_ID")
  # drop the site with the NA value (AQUE_2019-07-23) - just need to change the date
  # drop_na() %>% 


# check to make sure rows match up
# rownames(community_matrix) == community_metadata$sample_ID

# community matrix in presence-absence form
community_presabs <- community_matrix %>% 
  mutate(across(.cols = everything(), ~replace(., . > 0, 1)))

# 5. useful objects -------------------------------------------------------

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
                    "EGME", # Egregia menziesii
                    "Nandersoniana" # Nienburgia andersoniana
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



