# This is a script for loading in trait data.

############################################################-
# 0. source -------------------------------------------------
############################################################-

source(here::here("code", "00a_libraries.R"))

############################################################-
# 1. FvFm data ----------------------------------------------
############################################################-

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
  "20220824-CARP_MOHK-cleaned.csv",
  "20230620-NAPL-cleaned.csv",
  "20230627-MOHK-cleaned.csv",
  "20230627-MOHK-part2-cleaned.csv",
  "20230711-NAPL-cleaned.csv",
  "20230725-BULL-cleaned.csv",
  "20230731-CARP-cleaned.csv",
  "20230807-MOHK-cleaned.csv"))

# read in all .csv files into one data frame
fvfm_raw <- path %>% 
  map_df(~ read_csv(.))

############################################################-
# 2. other continuous traits --------------------------------
############################################################-

# generated csvs in `00b_google-sheets.R`
# remember to update file name if any have changed

# ⊣ a. metadata for subsamples ------------------------------

metadata_sub <- read_csv(here::here("data", "google-sheet-traits", "metadata_sub_sheet_2023-11-08.csv"))

# ⊣ b. metadata for individuals -----------------------------

metadata_ind <- read_csv(here::here("data", "google-sheet-traits", "metadata_ind_sheet_2023-11-08.csv"))

# ⊣ c. max height -------------------------------------------

ind_height <- read_csv(here::here("data", "google-sheet-traits", "ind_height_sheet_2023-11-08.csv"))

# ⊣ d. thallus length and width -----------------------------

lw <- read_csv(here::here("data", "google-sheet-traits", "lw_sheet_2023-11-08.csv")) 

# ⊣ e. thallus thickness ------------------------------------

thickness <- read_csv(here::here("data", "google-sheet-traits", "thickness_sheet_2023-11-08.csv"))

# ⊣ f. wet and dry weight -----------------------------------

weight <- read_csv(here::here("data", "google-sheet-traits", "weight_sheet_2024-09-26.csv"))
egme_weight <- read_csv(here::here("data", "google-sheet-traits", "egme_weight_sheet_2024-10-10.csv"))

# ⊣ g. surface area and perimeter ---------------------------

sa_peri <- read_csv(here::here("data", "google-sheet-traits", "sa_peri_sheet_2024-11-14.csv"))

# ⊣ h. branching order (still need to clean up) -------------

bra_ord <- read_csv(here::here("data", "google-sheet-traits", "bra_ord_sheet_2023-02-02.csv"))

# ⊣ i. toughness (still need to clean up) -------------------

toughness <- read_csv(here::here("data", "google-sheet-traits", "toughness_sheet_2023-02-02.csv")) 

# ⊣ j. volume -----------------------------------------------

volume <- read_csv(here::here("data", "google-sheet-traits", "volume_sheet_2024-06-20.csv"))

# ⊣ k. chlorophyll A ----------------------------------------

chlA <- read_csv(here::here("data", "google-sheet-traits", "chlA_sheet_2023-07-26.csv"))


# ⊣ l. isotopes ---------------------------------------------

isotopes <- readxl::read_xls(here("data", "isotopes", "Bui6695_Final.xls"), sheet = "clean",
                             na = "NA")

############################################################-
# 3. categorical traits -------------------------------------
############################################################-

coarse_traits <- read_csv(here::here("data", "algae-traits_literature-search_2022-02-28.csv"))

joe_traits <- read_csv(here::here("data", "fong-categorical", "Traits.csv")) %>% 
  mutate(species = str_to_sentence(species))

############################################################-
# 4. LTE data -----------------------------------------------
############################################################-

lte <- read_csv(here("data", "SBC-LTER-benthics", "LTE_All_Species_Biomass_at_transect_20220208.csv")) %>% 
  clean_names() %>% 
  mutate(across(c(site, treatment, group, mobility, growth_morph), tolower))

############################################################-
# 5. categorical trait data ---------------------------------
############################################################-

# in dropbox: /Users/An/Dropbox/grad-work/01_research/02-phd/03-functional-responses/data

cat_data <- repmis::source_data("https://www.dropbox.com/scl/fi/vpnue7o7dz2so6z1sxpwj/joe-traits-lter_2024-06-04.csv?rlkey=youehrjjnbunfc93gegiemopq&dl=1")




