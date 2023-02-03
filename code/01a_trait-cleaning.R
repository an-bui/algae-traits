# This is a script for cleaning up trait data.

############################################################-
# 0. source -------------------------------------------------
############################################################-

source(here::here("code", "00c_trait-data.R"))

############################################################-
# 1. cleaning/summarizing -----------------------------------
############################################################-

# ⊣ a. FvFm -------------------------------------------------

# subsample
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

# individual
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

# ⊣ b. thickness --------------------------------------------

# subsample
thickness_sub <- thickness %>% 
  pivot_longer(cols = thickness_01:thickness_10, names_to = "measurement_n", values_to = "thickness_mm") %>% 
  group_by(subsample_ID) %>% 
  summarize(thickness_mm_mean = mean(thickness_mm, na.rm = TRUE),
            thickness_mm_se = se(thickness_mm))

# individual
thickness_ind <- thickness %>% 
  pivot_longer(cols = thickness_01:thickness_10, names_to = "measurement_n", values_to = "thickness_mm") %>% 
  group_by(specimen_ID) %>% 
  summarize(thickness_mm_mean = mean(thickness_mm, na.rm = TRUE),
            thickness_mm_se = se(thickness_mm))

# ⊣ c. weights ----------------------------------------------

# subsample
weight_sub <- weight %>%
  select(specimen_ID, subsample_ID, weight_wet_mg, weight_wet_g, weight_dry_mg, weight_dry_g) %>%
  # calculate dry matter content
  mutate(dmc = weight_dry_mg/weight_wet_mg) %>%
  # take out specimen_ID column
  select(-specimen_ID)

# individual
weight_ind <- weight %>% 
  group_by(specimen_ID) %>% 
  summarize(total_wet = sum(weight_wet_mg, na.rm = TRUE),
            total_dry = sum(weight_dry_mg, na.rm = TRUE)) %>% 
  ungroup()

# ⊣ d. volume -----------------------------------------------

# subsample
volume_sub <- volume %>%
  select(subsample_ID, volume_total_mL)

# individual
volume_ind <- volume %>% 
  group_by(specimen_ID) %>% 
  summarize(total_volume = sum(volume_total_mL, na.rm = TRUE))

# ⊣ e. surface area and perimeter ---------------------------

# subsample
sa_peri_sub <- sa_peri %>% 
  select(specimen_ID, subsample_ID, area_total, peri_total) %>% 
  group_by(subsample_ID) %>% 
  summarize(area_total = sum(area_total),
            peri_total = sum(peri_total)) 

# ⊣ f. thallus length and width -----------------------------

# subsample
lw_sub <- lw %>% 
  select(-specimen_ID, -5)

# ⊣ g. branching order (not used) ----------------------------

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

# ⊣ h. toughness (not used) ----------------------------------

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

############################################################-
# 2. trait by sample matrices -------------------------------
############################################################-

# intermediate data frame for trait by species matrix
ct_prep <- metadata_sub %>% 
  select(specimen_ID, sp_code, lifestage) %>% 
  left_join(., coarse_traits, by = "sp_code") %>% 
  # select(-subsample_ID) %>% 
  select(specimen_ID, sp_code, scientific_name, 
         growth_form, pigment_type, life_habit, longevity, posture, branching_yn, lifestage) %>%
  unique()

# ⊣ a. leaf traits and average leaf values ------------------

# leaf traits
leaf_traits <- metadata_sub %>% 
  filter(type %in% c("whole", "thallus")) %>% 
  # filter(lifestage != "recruit") %>% 
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

# ⊣ b. individual values ------------------------------------

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

# ⊣ c. trait by species matrix ------------------------------

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

