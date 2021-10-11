
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
  # join with metadata
  left_join(., metadata, by = "specimen_ID") %>% 
  # select columns of interest
  select(specimen_ID, mean, var, sd, date_collected, site, sp_code) %>% 
  left_join(., algae_ct, by = "sp_code") %>% 
  unique()

# as a note: samples from 20210623 don't have FvFm, and one sample from IVEE doesn't either


# 3. thickness ------------------------------------------------------------

thickness_summary <- thickness %>% 
  pivot_longer(cols = thickness_01:thickness_10, names_to = "measurement_n", values_to = "thickness_mm") %>% 
  group_by(specimen_ID) %>% 
  summarize(mean = mean(thickness_mm),
            var = var(thickness_mm),
            sd = sd(thickness_mm)) %>% 
  # join with metadata
  left_join(., metadata_subsamples, by = "specimen_ID") %>% 
  # select columns of interest
  select(specimen_ID, mean, var, sd, date_collected, site, sp_code) %>% 
  left_join(., algae_ct, by = "sp_code") %>% 
  drop_na(sp_code)







