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

sa_peri <- read_csv(here::here("data", "google-sheet-traits", "sa_peri_sheet_2024-11-17.csv"))

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
# 4. LTER data ----------------------------------------------
############################################################-

lte <- read_csv(here("data", "SBC-LTER-benthics", "LTE_All_Species_Biomass_at_transect_20220208.csv")) %>% 
  clean_names() %>% 
  mutate(across(c(site, treatment, group, mobility, growth_morph), tolower))

biomass <- read_csv(here("data", "SBC-LTER-benthics",
                      "Annual_All_Species_Biomass_at_transect_20240823.csv")) %>% 
  clean_names() %>% 
  mutate(across(c(site, group, mobility, growth_morph), tolower)) %>% 
  filter(group == "algae") %>% 
  mutate(dry_gm2 = replace(dry_gm2, dry_gm2 == -99999, NA)) 

lter_algae_group_biomass <- biomass %>% 
  filter(sp_code != "MAPY") %>% 
  filter(!(site %in% c("scdi", "sctw"))) %>% 
  # filter(sp_code %in% algae_spcode_factors) %>% 
  left_join(., enframe(groups), by = c("sp_code" = "name")) %>% 
  replace_na(list(value = "other")) %>% 
  rename("new_groups" = "value") %>% 
  unite(col = "ID", date, site, transect, sep = "_", remove = FALSE) %>% 
  # filter(ID %in% fish_meta$ID) %>% 
  select(ID, year, month, date, site, transect, dry_gm2, new_groups) %>%  
  group_by(ID, new_groups) %>% 
  mutate(group_biomass = sum(dry_gm2, na.rm = TRUE)) %>% 
  ungroup() %>% 
  select(-dry_gm2) %>% 
  unique() %>% 
  group_by(ID) %>% 
  mutate(total_biomass = sum(group_biomass, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(proportion = group_biomass/total_biomass) %>% 
  select(!group_biomass) %>% 
  pivot_wider(names_from = new_groups, values_from = proportion) %>% 
  unique() %>% 
  left_join(., select(fish_richness, ID, fish_sppnum), by = "ID") %>% 
  clean_names()

lter_algae_group_biomass_filtered <- lter_algae_group_biomass %>% 
  select(id, foliose_bushy_turf:bushy_short_red) %>% 
  column_to_rownames("id")

algae_matrix <- biomass %>% 
  unite("ID", date, site, transect, remove = FALSE) %>% 
  select(ID, sp_code, dry_gm2) %>% 
  filter(sp_code %in% algae_spcode_factors) %>% 
  pivot_wider(names_from = sp_code, values_from = dry_gm2) %>% 
  column_to_rownames("ID") %>% 
  mutate(sum = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>% 
  filter(sum > 0) %>% 
  select(!sum)

group_df <- biomass %>% 
  select(sp_code) %>% 
  unique() %>% 
  filter(sp_code %in% algae_spcode_factors) %>% 
  left_join(., enframe(groups), by = c("sp_code" = "name")) %>% 
  replace_na(list(value = "other")) %>% 
  column_to_rownames("sp_code")

comm_fun <- FD::dbFD(x = group_df,
                     a = algae_matrix)

evenness <- enframe(comm_fun$FEve) %>% 
  rename("id" = "name",
         "even" = "value") %>% 
  left_join(., lter_algae_group_biomass, by = "id")

ggplot(data = evenness,
       aes(x = even,
           y = fish_sppnum)) +
  geom_point(aes(color = articulated_coralline,
                 size = articulated_coralline),
             shape = 21) +
  geom_smooth()

ggplot(data = lter_algae_group_biomass,
       aes(x = tall_stipate,
           y = fish_sppnum)) +
  geom_point() +
  geom_smooth(se = TRUE, method = "lm")

mod <- glmmTMB(fish_sppnum ~ tall_stipate + (1|year) + (1|site) + (1|transect),
           data = lter_algae_group_biomass,
           family = nbinom1(link = "log"))

plot(simulateResiduals(mod))
summary(mod)
plot(ggpredict(mod), show_data = TRUE)

irr <- read_csv(here("data", "SBC-LTER-benthics", "Hourly_Irrandiance_All_Year_20220131.csv")) %>% 
  clean_names() %>% 
  mutate(across(where(is.character), tolower))

fish <- read_csv(here("data", "SBC-LTER-benthics",
                      "Annual_fish_comb_20240823.csv")) %>% 
  clean_names() %>% 
  mutate(across(c(site, group, mobility, growth_morph), tolower)) %>% 
  unite(col = "ID", date, site, transect, sep = "_", remove = FALSE) %>% 
  filter(!(site %in% c("scdi", "sctw"))) %>% 
  mutate(count = replace(count, count == -99999, NA)) 

fish_matrix <- fish %>% 
  select(ID, scientific_name, count) %>% 
  group_by(ID, scientific_name) %>% 
  summarize(count = sum(count, na.rm = TRUE)) %>% 
  ungroup() %>% 
  complete(ID, scientific_name, fill = list(count = 0)) %>% 
  pivot_wider(names_from = scientific_name, values_from = count) %>% 
  mutate(sum = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>% 
  filter(sum > 0) %>% 
  select(-sum) %>% 
  column_to_rownames("ID")

fish_meta <- fish %>% 
  select(ID, year, month, date, site, transect) %>% 
  unique()

# 2001_8_bull_1 done twice? 2001-08-29 and 2001-08-17

setdiff(fish_meta$ID, lter_filtered_algae$ID)

fish_richness <- vegan::specnumber(fish_matrix) %>% 
  as_tibble(rownames = "ID") %>% 
  left_join(., fish_meta, by = "ID") %>% 
  rename("fish_sppnum" = "value")

ggplot(data = fish_richness %>% filter(year > 2020),
       aes(x = date,
           y = fish_sppnum,
           color = site)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~site)

zero_surveys <- fish %>% 
  select(ID, scientific_name, count) %>% 
  group_by(ID, scientific_name) %>% 
  summarize(count = sum(count, na.rm = TRUE)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = scientific_name, values_from = count) %>% 
  mutate(sum = rowSums(across(where(is.numeric)), na.rm = TRUE)) %>% 
  filter(sum == 0) %>% 
  pull(ID)

# irr %>% 
#   filter(site == "aque",
#          transect %in% c(4, 8)) %>% 
#   mutate(transect = case_when(
#     transect == 4 ~ "reference",
#     transect == 8 ~ "removal"
#   )) %>% 
#   group_by(sensor_location,
#            date_local, 
#            transect) %>% 
#   summarize(mean_light = mean(light_umol)) %>% 
#   ggplot(aes(x = date_local,
#              y = mean_light,
#              color = transect)) +
#   geom_point()

############################################################-
# 5. categorical trait data ---------------------------------
############################################################-

# in dropbox: /Users/An/Dropbox/grad-work/01_research/02-phd/03-functional-responses/data

cat_data <- repmis::source_data("https://www.dropbox.com/scl/fi/vpnue7o7dz2so6z1sxpwj/joe-traits-lter_2024-06-04.csv?rlkey=youehrjjnbunfc93gegiemopq&dl=1")




