# This is a script to get the range of weights for stable isotope samples.

############################################################-
# 0. source -------------------------------------------------
############################################################-

source(here::here("code", "00c_trait-data.R"))

############################################################-
# 1. data ---------------------------------------------------
############################################################-

biomass_rels <- read_csv(here("data", "isotopes", "Algal_Biomass_Relationships_NPP_calculation_20210113.csv")) %>% 
  clean_names() %>% 
  mutate(sp_code = case_when(
    scientific_name == "Nienburgia andersoniana" ~ "Nandersoniana",
    TRUE ~ sp_code
  ))

############################################################-
# 2. calculating sample weights -----------------------------
############################################################-

# ⊣ a. sample weights ---------------------------------------

si_weights <- biomass_rels %>% 
  # calculate low and high ranges for delta C and N
  mutate(c_sample_weight_low_mg = (0.071*100)/carbon,
         c_sample_weight_high_mg = (0.64*100)/carbon,
         n_sample_weight_low_mg = (0.01*100)/nitrogen,
         n_sample_weight_high_mg = (0.09*100)/nitrogen,
         # calculate ideal sample weights based on nitrogen ranges
         ideal_sample_weight_low_mg = (n_sample_weight_low_mg + n_sample_weight_high_mg)/2,
         ideal_sample_weight_high_mg = ideal_sample_weight_low_mg + ideal_sample_weight_low_mg*0.2) %>% 
  # round to 2 decimal points
  mutate(across(c_sample_weight_low_mg:ideal_sample_weight_high_mg, ~ round(.x, 2)))

# ⊣ b. sample lists -----------------------------------------

si_samples_run1 <- c(
  "20230620-NAPL-003-B",
  "20230620-NAPL-005-D",
  "20230620-NAPL-010-B",
  "20230620-NAPL-015-B",
  "20230620-NAPL-016-A",
  "20230620-NAPL-022-A",
  "20230620-NAPL-029-A",
  "20230620-NAPL-032-A",
  "20230627-MOHK-001-C",
  "20230627-MOHK-004-C",
  "20230627-MOHK-012-D",
  "20230627-MOHK-013-B",
  "20230627-MOHK-016-B",
  "20230627-MOHK-020-D",
  "20230627-MOHK-024-A",
  "20230711-NAPL-006-A",
  "20230711-NAPL-010-B",
  "20230711-NAPL-014-A",
  "20230711-NAPL-015-B",
  "20230711-NAPL-018-A",
  "20230711-NAPL-021-D",
  "20230725-BULL-002-A",
  "20230725-BULL-004-C",
  "20230725-BULL-009-B",
  "20230725-BULL-011-C",
  "20230731-CARP-003-D",
  "20230731-CARP-005-B",
  "20230731-CARP-010-A",
  "20230731-CARP-012-B",
  "20230807-MOHK-001-C",
  "20230807-MOHK-010-A",
  "20230807-MOHK-012-A",
  "20230807-MOHK-015-C",
  "20230807-MOHK-020-A",
  "20230807-MOHK-025-A"
  )

si_samples <- metadata_sub %>% 
  dplyr::filter(subsample_ID %in% si_samples_run1) %>% 
  left_join(., si_weights, by = "sp_code") %>% 
  select(subsample_ID, sp_code, c_sample_weight_low_mg:ideal_sample_weight_high_mg)

write_csv(si_samples, file = here("data", "isotopes", "si_run01.csv"))

