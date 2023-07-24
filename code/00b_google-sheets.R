# This is a script to get data from google drive and save each data frame as a csv. 
# Note to self: THIS ONLY NEEDS TO BE RUN IF YOU HAVE UPDATED ANY OF THE DATA.

############################################################-
# 1. function to write csv ----------------------------------
############################################################-

sheet_write_csv <- function(df) {
  # convert object name to string
  df_name <- deparse(substitute(df))
  # write csv from object and put into google-sheet-traits folder
  write_csv(df, here::here("data", "google-sheet-traits", paste(df_name, "_", today(), ".csv", sep = "")))
}

############################################################-
# 2. data from google sheet ---------------------------------
############################################################-

# packages
source(here::here("code", "00a_libraries.R"))

# file id from url
sheet_id <- "1h2eHoL5kXMRwExt3hjrLavxG0pajHOvA1UfJd24pRk4"

# googlesheets4 will ask for authentication from whichever read_sheet command runs first

# ⊣ a. metadata for subsamples ------------------------------

metadata_sub_sheet <- read_sheet(sheet_id, sheet = "00b-metadata_sub") %>% 
  mutate(date_collected = ymd(date_collected)) %>% 
  mutate(year = year(date_collected))

# last updated: 2023-07-24
sheet_write_csv(metadata_sub_sheet)

# ⊣ b. metadata for individuals -----------------------------

metadata_ind_sheet <- metadata_sub_sheet %>% 
  select(specimen_ID, date_collected, year, site, sp_code, lifestage) %>% 
  unique()

# last updated: 2023-07-24
sheet_write_csv(metadata_ind_sheet)

# ⊣ c. max height -------------------------------------------

ind_height_sheet <- read_sheet(sheet_id, sheet = "02a-ind_height", na = "NA") 

# last updated: 2023-02-02
sheet_write_csv(ind_height_sheet)

# ⊣ d. thallus length and width -----------------------------

lw_sheet <- read_sheet(sheet_id, sheet = "02b-lw", na = "NA") 

# last updated: 2023-02-02
sheet_write_csv(lw_sheet)

# ⊣ e. thallus thickness ------------------------------------

thickness_sheet <- read_sheet(sheet_id, sheet = "03-thickness", na = "NA")

# last updated: 2023-02-02
sheet_write_csv(thickness_sheet)

# ⊣ f. wet and dry weight -----------------------------------

weight_sheet <- read_sheet(sheet_id, sheet = "04-weight", na = "NA")

# last updated: 2023-02-02
sheet_write_csv(weight_sheet)

# ⊣ g. surface area and perimeter ---------------------------

sa_peri_sheet <- read_sheet(sheet_id, sheet = "05a-scans") 

# last updated: 2023-02-20
sheet_write_csv(sa_peri_sheet)

# ⊣ h. branching order (still need to clean up) -------------

bra_ord_sheet <- read_sheet(sheet_id, sheet = "05b-branching_order", na = "NA") 

# last updated: 2023-02-02
sheet_write_csv(bra_ord_sheet)

# ⊣ i. toughness (still need to clean up) -------------------

toughness_sheet <- read_sheet(sheet_id, sheet = "06-toughness", na = "NA") 

# last updated: 2023-02-02
sheet_write_csv(toughness_sheet)

# ⊣ j. volume -----------------------------------------------

volume_sheet <- read_sheet(sheet_id, sheet = "07-volume", na = "NA") 

# last updated: 2023-04-06
sheet_write_csv(volume_sheet)

# ⊣ k. chlorophyll A -----------------------------------------

chlA_sheet <- read_sheet(sheet_id, sheet = "08-chlA", na = "NA") 

# last updated: 2023-07-23
sheet_write_csv(chlA_sheet)



