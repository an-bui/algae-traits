
# 1.  libraries -----------------------------------------------------------

library(tidyverse)
library(googlesheets4)
library(here)

# 2.  getting data from google drive --------------------------------------

# gets the file data from google drive, only have to do this once to get the sheet id
# trait_data_id <- googledrive::drive_get("trait_data")

sheet_id <- "1h2eHoL5kXMRwExt3hjrLavxG0pajHOvA1UfJd24pRk4"

metadata <- read_sheet(sheet_id, sheet = "00b-metadata")

hw <- read_sheet(sheet_id, sheet = "02-hw")

thickness <- read_sheet(sheet_id, sheet = "03-thickness")

weight <- read_sheet(sheet_id, sheet = "04-weight")

sa_peri <- read_sheet(sheet_id, sheet = "05a-scans")

bra_ord <- read_sheet(sheet_id, sheet = "05b-branching_order")

toughness <- read_sheet(sheet_id, sheet = "06-toughness")

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



