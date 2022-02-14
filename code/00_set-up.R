
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

## a. traits 
# coarse trait data
coarse_traits <- read_csv(here::here("data", "algae_all.csv"))

## b. benthics
benthic_cleaning_fxn <- function(df) {
  df %>% 
    clean_names() %>% 
    # create a sample_ID for each sampling date at each site
    unite("sample_ID", site, date, remove = FALSE) %>% 
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

# biomass
biomass <- read_csv(here::here("data", "SBC-LTER-benthics", 
                               "Annual_All_Species_Biomass_at_transect_20210108.csv")) %>% 
  benthic_cleaning_fxn() %>% 
  # replace all -99999 values with 0
  mutate(dry_gm2 = replace(dry_gm2, dry_gm2 < 0, 0)) %>% 
  mutate(date = ymd(date))

# percent cover
percov <- read_csv(here::here("data", "SBC-LTER-benthics", 
                              "Annual_Cover_All_Years_20210108.csv")) %>% 
  benthic_cleaning_fxn()


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



