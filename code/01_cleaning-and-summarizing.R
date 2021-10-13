
# 1. set up ---------------------------------------------------------------

source(here::here("code", "00_set-up.R"))


# 2. FvFm -----------------------------------------------------------------

fvfm_summary <- fvfm_raw %>% 
  # only use FvFm, drop the NAs, change the column name
  select(specimen_ID, '1:Fv/Fm') %>% 
  drop_na() %>%
  rename('fvfm_meas' = '1:Fv/Fm') %>% 
  # make sure that fvfm_meas is numeric
  mutate(fvfm_meas = as.numeric(fvfm_meas)) %>%
  # calculate mean, variance, standard deviation
  group_by(specimen_ID) %>% 
  summarize(mean = mean(fvfm_meas),
          var = var(fvfm_meas),
          sd = sd(fvfm_meas)) %>% 
  # join with metadata: spp codes + specimen_ID
  left_join(., metadata, by = "specimen_ID") %>% 
  # select columns of interest
  select(specimen_ID, mean, var, sd, date_collected, site, sp_code) %>% 
  # join with algae_ct: spp codes + taxonomic info
  left_join(., algae_ct, by = "sp_code") %>% 
  unique()

# as a note: samples from 20210623 don't have FvFm, and one sample from IVEE doesn't either


# 3. thickness ------------------------------------------------------------

# 20210709-MOHK-004 only has 7 measurements for thickness, for whatever reason
# this calculates the average thickness of those 7 measurements
quickfix <- thickness %>% 
  filter(specimen_ID == "20210709-MOHK-004") %>% 
  select(1:8) %>% 
  pivot_longer(cols = 2:8, names_to = "measurement_n", values_to = "thickness_mm") %>% 
  pull(thickness_mm) %>% 
  mean()

# some problems with the subsampling - used the metadata_subsamples df

thickness_summary <- thickness %>% 
  pivot_longer(cols = thickness_01:thickness_10, names_to = "measurement_n", values_to = "thickness_mm") %>% 
  # replace NAs for 20210709-MOHK-004 with the quickfix value
  mutate(thickness_mm = replace(thickness_mm, 
                                is.na(thickness_mm) & specimen_ID == "20210709-MOHK-004", 
                                quickfix)) %>% 
  group_by(specimen_ID) %>% 
  summarize(mean = mean(thickness_mm),
            var = var(thickness_mm),
            sd = sd(thickness_mm)) %>% 
  # join with metadata_subsamples: spp codes + specimen_ID
  left_join(., metadata_subsamples, by = "specimen_ID") %>% 
  # select columns of interest
  select(specimen_ID, mean, var, sd, date_collected, site, sp_code) %>% 
  # join with algae_ct: spp codes + taxonomic info
  left_join(., algae_ct, by = "sp_code") %>% 
  drop_na(sp_code)

# 4. dry weight -----------------------------------------------------------

# only whole dry weights

dryw_summary <- weight %>% 
  filter(type == "whole") %>% 
  # join with metadata: spp codes + specimen_ID
  left_join(., metadata, by = "specimen_ID") %>% 
  drop_na(sp_code) %>% 
  # join with algae_ct: spp codes + taxonomic info
  left_join(., algae_ct, by = "sp_code")



# 5. volume ---------------------------------------------------------------

volume_summary <- volume %>% 
  filter(type == "whole") %>% 
  # join with metadata: spp codes + specimen_ID
  left_join(., metadata_subsamples, by = "specimen_ID") %>% 
  drop_na(sp_code) %>% 
  # join with algae_ct: spp codes + taxonomic info
  left_join(., algae_ct, by = "sp_code")


# 6. surface area ---------------------------------------------------------

sa_summary <- sa_peri %>% 
  # drop all observations that haven't been processed
  drop_na(results_file) %>% 
  select(specimen_ID, area_total) %>% 
  group_by(specimen_ID) %>% 
  summarize(area_total = sum(area_total)) %>% 
  drop_na(area_total) %>% 
  # join with metadata_subsamples: spp code + specimen_ID
  left_join(., metadata_subsamples, by = "specimen_ID") %>% 
  # select columns of interest
  select(specimen_ID, area_total, date_collected, site, sp_code) %>% 
  # join with algae_ct: spp codes + taxonomic info
  left_join(., algae_ct, by = "sp_code")


# 7. perimeter ------------------------------------------------------------

peri_summary <- sa_peri %>% 
  # drop all observations that haven't been processed
  drop_na(results_file) %>% 
  select(specimen_ID, peri_total) %>% 
  group_by(specimen_ID) %>% 
  summarize(peri_total = sum(peri_total)) %>% 
  # join with metadata_subsamples: spp code + specimen_ID
  left_join(., metadata_subsamples, by = "specimen_ID") %>% 
  # select columns of interest
  select(specimen_ID, peri_total, date_collected, site, sp_code) %>% 
  # join with algae_ct: spp codes + taxonomic info
  left_join(., algae_ct, by = "sp_code")


# 8. surface area:perimeter ratio -----------------------------------------

ratio_saperi <- sa_summary %>% 
  select(specimen_ID, area_total) %>% 
  full_join(., peri_summary, by = "specimen_ID") %>% 
  mutate(ratio = area_total/peri_total) %>% 
  relocate(specimen_ID, area_total, peri_total, ratio)


# 9. surface area:volume ratio --------------------------------------------

ratio_savolume <- sa_summary %>% 
  select(specimen_ID, area_total) %>% 
  full_join(., volume_summary, by = "specimen_ID") %>% 
  # drop specimens that don't have area
  drop_na(area_total) %>% 
  mutate(ratio = area_total/volume_total_mL) %>% 
  drop_na(ratio)

# 10. max height and width ------------------------------------------------

hw_summary <- hw %>% 
  filter(type == "whole") %>% 
  # join with metadata: spp codes + specimen_ID
  left_join(., metadata_subsamples, by = "specimen_ID") %>% 
  drop_na(sp_code) %>% 
  # join with algae_ct: spp codes + taxonomic info
  left_join(., algae_ct, by = "sp_code")

