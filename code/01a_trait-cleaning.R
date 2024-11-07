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
  filter(fvfm_meas != "-") %>% 
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

holdfasts <- metadata_sub %>% 
  filter(type %in% c("holdfast")) %>% 
  pull(subsample_ID)

stipes_holdfasts <- metadata_sub %>% 
  filter(type %in% c("stipe", "holdfast")) %>% 
  pull(subsample_ID)

# subsample
thickness_sub <- thickness %>% 
  # take out stipes and holdfasts
  filter(!(subsample_ID %in% c(holdfasts))) %>% 
  pivot_longer(cols = thickness_01:thickness_10, names_to = "measurement_n", values_to = "thickness_mm") %>% 
  group_by(subsample_ID) %>% 
  summarize(thickness_mm_mean = mean(thickness_mm, na.rm = TRUE),
            thickness_mm_se = se(thickness_mm))

# individual
thickness_ind <- thickness %>% 
  # take out stipes and holdfasts
  filter(!(subsample_ID %in% c(holdfasts))) %>% 
  pivot_longer(cols = thickness_01:thickness_10, names_to = "measurement_n", values_to = "thickness_mm") %>% 
  group_by(specimen_ID) %>% 
  summarize(thickness_mm_mean = mean(thickness_mm, na.rm = TRUE),
            thickness_mm_se = se(thickness_mm))

# ⊣ c. weights ----------------------------------------------

# subsample
weight_sub <- weight %>%
  # take out stipes and holdfasts
  filter(!(subsample_ID %in% c(stipes_holdfasts))) %>% 
  select(specimen_ID, subsample_ID, weight_wet_mg, weight_wet_g, weight_dry_mg, weight_dry_g) %>%
  # calculate dry matter content
  mutate(dmc = weight_dry_mg/weight_wet_mg) %>%
  # take out specimen_ID column
  select(-specimen_ID)

# EGME weights
egme_weight

# individual
weight_ind <- weight %>% 
  group_by(specimen_ID) %>% 
  summarize(total_wet = sum(weight_wet_mg, na.rm = TRUE),
            total_dry = sum(weight_dry_mg, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(total_dmc = total_dry/total_wet) %>% 
  mutate(total_dry = case_when(
    specimen_ID == pluck(egme_weight, 1, 1) ~ pluck(egme_weight, 4, 1),
    specimen_ID == pluck(egme_weight, 1, 2) ~ pluck(egme_weight, 4, 2),
    specimen_ID == pluck(egme_weight, 1, 3) ~ pluck(egme_weight, 4, 3),
    specimen_ID == pluck(egme_weight, 1, 4) ~ pluck(egme_weight, 4, 4),
    specimen_ID == pluck(egme_weight, 1, 5) ~ pluck(egme_weight, 4, 5),
    TRUE ~ total_dry
  ))

# ⊣ d. volume -----------------------------------------------

# subsample
volume_sub <- volume %>%
  select(subsample_ID, volume_total_mL) %>% 
  # for samples that had too small of a volume to read, make them have 1 mL volume
  mutate(volume_total_mL = case_when(
    subsample_ID %in% c("20210719-IVEE-009-A", "20210719-IVEE-010-A", "20220407-AQUE-010-F") ~  1,
    TRUE ~ volume_total_mL
  ))

# individual
volume_ind <- volume %>% 
  # for samples that had too small of a volume to read, make them have 1 mL volume
  mutate(volume_total_mL = case_when(
    subsample_ID %in% c("20210719-IVEE-009-A", "20210719-IVEE-010-A", "20220407-AQUE-010-F") ~  1,
    TRUE ~ volume_total_mL
  )) %>% 
  group_by(specimen_ID) %>% 
  summarize(total_volume = sum(volume_total_mL, na.rm = TRUE))

# ⊣ e. surface area and perimeter ---------------------------

pneumatocysts <- sa_peri %>% 
  filter(notes_scans == "pneumatocyst") %>% 
  pull(subsample_ID)

# subsample
sa_peri_sub <- sa_peri %>% 
  # take out measurements for EGME pneumatocysts
  filter(!(subsample_ID %in% pneumatocysts)) %>% 
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

# ⊣ i. chlorophyll A ----------------------------------------

chlA_sub <- chlA %>% 
  # calculating total chlA ug per mg wet weight
  mutate(chlA_ug = 
           # chlorophyll in DMSO
           (2 * dmso_ug_L/1000 * (dmso_sample + dmso_diluent)/dmso_sample) +
           # chlorophyll in MilliQ wash
           (2 * milliq_ug_L/1000) +
           # chlorophyll in 3:1:1 acetone:methanol:MilliQ
           (5 * ace_met_mil_ug_L/1000 * (ace_met_mil_sample + ace_met_mil_diluent)/ace_met_mil_sample),
         chlA_ug_mg2 = chlA_ug/wet_weight_mg
           )

chlA_ind <- chlA_sub %>% 
  group_by(specimen_ID) %>% 
  summarize(chlA_mean = mean(chlA_ug_mg, na.rm = TRUE),
            chlA_se = se(chlA_ug_mg)) %>% 
  ungroup()

# ⊣ j. isotopes ---------------------------------------------

isotopes_sub <- isotopes %>% 
  filter(!is.na(subsample_ID)) %>% 
  select(!c(sample_id_run, amount_ug)) 

isotopes_ind <- isotopes_sub %>% 
  separate_wider_delim(cols = subsample_ID, names = c("date", "site", "ind", "sub"), delim = "-", cols_remove = FALSE) %>% 
  unite("specimen_ID", date:ind, sep = "-") %>% 
  select(!sub) %>% 
  group_by(specimen_ID) %>% 
  summarize(delta_15_n_mean = mean(delta_15_n, na.rm = TRUE),
            delta_15_n_se = se(delta_15_n),
            delta_13_c_mean = mean(delta_13_c, na.rm = TRUE),
            delta_13_c_se = se(delta_13_c),
            weight_percent_n_mean = mean(weight_percent_n, na.rm = TRUE),
            weight_percent_n_se = se(weight_percent_n),
            weight_percent_c_mean = mean(weight_percent_c, na.rm = TRUE),
            weight_percent_c_se = se(weight_percent_c))

############################################################-
# 2. trait by sample matrices -------------------------------
############################################################-

# vector of recruit subsamples
recruits <- metadata_sub %>% 
  filter(lifestage == "recruit") %>% 
  pull(specimen_ID)

# intermediate data frame for trait by species matrix
ct_prep <- metadata_sub %>% 
  # filter(!(specimen_ID %in% recruits)) %>% 
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
  filter(!(specimen_ID %in% recruits)) %>% 
  left_join(., fvfm_sub, by = "subsample_ID") %>% 
  left_join(., thickness_sub, by = "subsample_ID") %>% 
  left_join(., weight_sub, by = "subsample_ID") %>% 
  left_join(., volume_sub, by = "subsample_ID") %>% 
  left_join(., lw_sub, by = "subsample_ID") %>% 
  left_join(., sa_peri_sub, by = "subsample_ID") %>% 
  left_join(., isotopes_sub, by = "subsample_ID") %>% 
  # ratios
  mutate(sap_ratio = area_total/peri_total,
         sav_ratio = area_total/volume_total_mL,
         aspect_ratio = length/width,
         sta_mm_mg = area_total/weight_dry_mg)

# average leaf values
av_leaf_values <- leaf_traits %>% 
  group_by(specimen_ID) %>% 
  summarize(sta_mean = mean(sta_mm_mg, na.rm = TRUE),
            sta_se = se(sta_mm_mg),
            frond_dmc_mean = mean(dmc, na.rm = TRUE),
            frond_dmc_se = se(dmc),
            sav_mean = mean(sav_ratio, na.rm = TRUE),
            sav_se = se(sav_ratio),
            sap_mean = mean(sap_ratio, na.rm = TRUE),
            sap_se = se(sap_ratio),
            aspect_ratio_mean = mean(aspect_ratio, na.rm = TRUE),
            aspect_ratio_se = se(aspect_ratio),
            frond_length_mean = mean(length, na.rm = TRUE),
            frond_length_se = se(length),
            frond_width_mean = mean(width, na.rm = TRUE), 
            frond_width_se = se(width),
            frond_area_mean = mean(area_total, na.rm = TRUE),
            frond_area_se = se(area_total),
            frond_peri_mean = mean(peri_total, na.rm = TRUE),
            frond_peri_se = se(peri_total),
            frond_ww_mean = mean(weight_wet_mg, na.rm = TRUE),
            frond_ww_se = se(weight_wet_mg),
            frond_dw_mean = mean(weight_dry_mg, na.rm = TRUE),
            frond_dw_se = se(weight_dry_mg),
            frond_volume_mean = mean(volume_total_mL, na.rm = TRUE),
            frond_volume_se = se(volume_total_mL))

# ⊣ b. individual values ------------------------------------

# all traits for each individual with mean taken for each trait
ind_traits <- ct_prep %>% 
  left_join(., ind_height, by = "specimen_ID") %>% 
  left_join(., thickness_ind, by = "specimen_ID") %>% 
  left_join(., av_leaf_values, by = "specimen_ID") %>% 
  left_join(., fvfm_ind, by = "specimen_ID") %>% 
  left_join(., weight_ind, by = "specimen_ID") %>% 
  left_join(., volume_ind, by = "specimen_ID") %>% 
  left_join(., chlA_ind, by = "specimen_ID") %>% 
  left_join(., isotopes_ind, by = "specimen_ID") %>% 
  mutate(mass_to_height = total_dry/maximum_height) %>% 
  left_join(., (metadata_ind %>% select(specimen_ID, date_collected, site)), by = "specimen_ID") %>% 
  filter(sp_code != "EGME") %>% 
  filter(!(sp_code == "PTCA" & lifestage == "recruit")) %>% 
  filter(!(specimen_ID %in% c("20220726-BULL-016",
                              "20220726-BULL-024",
                              "20220421-AQUE-005"))) %>% 
  mutate(date_collected = ymd(date_collected)) %>% 
  mutate(year = year(date_collected)) %>% 
  mutate(sav_scaled = sav_mean*total_volume,
         frond_area_scaled = frond_area_mean*total_wet,
         sta_scaled = sta_mean*total_dry,
         frond_peri_scaled = frond_peri_mean*total_wet,
         frond_length_scaled = frond_length_mean*total_wet,
         frond_width_scaled = frond_width_mean*total_wet)

# ⊣ c. trait by species matrix ------------------------------

# trait values averaged across sites
# tbspp_matrix <- metadata_ind %>% 
#   select(specimen_ID) %>% 
#   unique() %>% 
#   left_join(., ind_height, by = "specimen_ID") %>% 
#   left_join(., thickness_ind, by = "specimen_ID") %>% 
#   left_join(., av_leaf_values, by = "specimen_ID") %>% 
#   left_join(., fvfm_ind, by = "specimen_ID") %>% 
#   left_join(., weight_ind, by = "specimen_ID") %>% 
#   left_join(., volume_ind, by = "specimen_ID") %>% 
#   left_join(., metadata_ind, by = "specimen_ID") %>% 
#   filter(!(sp_code == "PTCA" & lifestage == "recruit")) %>% 
#   select(specimen_ID, sp_code,
#          fvfm_mean, thickness_mm_mean, total_dry, tdmc_mean,
#          total_volume, sta_mean, sap_mean, sav_mean, maximum_height) %>% 
#   # calculate mean trait value for each species
#   group_by(sp_code) %>% 
#   summarize_at(vars(fvfm_mean:maximum_height), mean, na.rm = TRUE) %>% 
#   ungroup() %>% 
#   drop_na(sp_code) %>% 
#   # join with coarse traits
#   left_join(., coarse_traits, by = "sp_code") %>% 
#   select(sp_code:taxon_phylum, growth_form, pigment_type, life_habit, longevity, posture, branching_yn) %>% 
#   filter(sp_code != "MAPY") %>% 
#   filter(sp_code %in% algae_proposal) %>% 
#   column_to_rownames("sp_code") %>% 
#   # replace NaNs with NAs
#   mutate_all(na_if, "NaN")




############################################################-
# 3. categorical JoE traits ---------------------------------
############################################################-

lte_spp <- lte %>% 
  filter(group == "algae") %>% 
  filter(!scientific_name %in% c("Unidentifiable juvenile kelp", 
                              # "Halymenia spp.; Schizymenia pacifica",
                              "Unidentifiable Branching Red Alga", 
                              "small Ceramiaceae spp.",
                              "Unidentifiable small brown blade")) %>% 
  select(scientific_name, taxon_phylum, taxon_order, taxon_family) %>% 
  unique() %>% 
  left_join(., joe_traits, by = c("scientific_name" = "species"))
# 14 spp have traits already from Fong et al. JoE, 40 species do not, 1 "species" combines two genera

# write_csv(lte_spp, file = here("data", "fong-categorical", "joe-traits-lter.csv"))













