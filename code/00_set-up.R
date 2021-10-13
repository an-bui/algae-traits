
# 1.  libraries -----------------------------------------------------------

library(tidyverse)
library(googlesheets4)
library(here)

# 2.  getting data from google drive --------------------------------------

# gets the file data from google drive
trait_data_id <- googledrive::drive_get("trait_data")

sheet_id <- "1h2eHoL5kXMRwExt3hjrLavxG0pajHOvA1UfJd24pRk4"

metadata <- read_sheet(sheet_id, sheet = "00b-metadata") %>% 
  filter(site != "Campus Point")

metadata_subsamples <- read_sheet(sheet_id, sheet = "00c-metadata-subsamples") %>% 
  filter(site != "Campus Point")

hw <- read_sheet(sheet_id, sheet = "02-hw") 

thickness <- read_sheet(sheet_id, sheet = "03-thickness", na = "NA") 

weight <- read_sheet(sheet_id, sheet = "04-weight", na = "NA") 

sa_peri <- read_sheet(sheet_id, sheet = "05a-scans") 

bra_ord <- read_sheet(sheet_id, sheet = "05b-branching_order") 

toughness <- read_sheet(sheet_id, sheet = "06-toughness", na = "NA") 

volume <- read_sheet(sheet_id, sheet = "07-volume") 


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

# coarse trait data
coarse_traits <- read_csv(here::here("data", "00-coarse_traits.csv"))

# all species with group, mobility, species code, scientific name
algae_all <- read_csv(here::here("data", "spp_names.csv")) %>% 
  filter(group == "algae") %>% 
  select(-group_mobility)


# 5.  useful vectors and data frames --------------------------------------

# most abundant algae
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

# data frame with scientific name, taxonomy, and coarse traits
# algae_ct = "algae coarse traits"
algae_ct <- full_join(algae_all, coarse_traits, by = "sp_code") %>% 
  drop_na(sp_code)


