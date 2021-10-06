
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
  full_join(., metadata, by = "specimen_ID") %>% 
  # select columns of interest
  select(specimen_ID, mean, var, sd, date_collected, site, species)
