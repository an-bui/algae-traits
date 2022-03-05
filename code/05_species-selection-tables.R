table_df <- biomass %>% 
  select(site, scientific_name, dry_gm2, wm_gm2) %>% 
  filter(site %in% c("napl")) %>% 
  # mutate(site = fct_relevel(site, "aque", "napl", "ivee", "mohk", "carp")) %>% 
  # filter out MAPY
  filter(scientific_name != "Macrocystis pyrifera") %>% 
  # create a new column where total biomass across sites for the whole survey is calculated
  # group_by(site) %>% 
  mutate(total_dry = sum(dry_gm2),
         total_wet = sum(wm_gm2)) %>% 
  # ungroup() %>% 
  # create a new column where total biomass for the species for the whole survey is calculated
  group_by(scientific_name) %>% 
  summarize(total_sp_dry = sum(dry_gm2),
         total_sp_wet = sum(wm_gm2), 
         percent_sp_dry = total_sp_dry/total_dry,
         percent_sp_wet = total_sp_wet/total_wet) %>% 
  ungroup() %>% 
  unique() %>% 
  select(scientific_name, percent_sp_dry) %>% 
  mutate(percent_sp_dry = round(percent_sp_dry, 4)) %>% 
  mutate(percent_whole_sp_dry = percent_sp_dry*100) %>% 
  select(-percent_sp_dry) %>% 
  left_join(., coarse_traits, by = "scientific_name") 
  filter(taxon_phylum == "Rhodophyta")
  select(scientific_name, sp_code, site, percent_whole_sp_dry) %>% 
  pivot_wider(names_from = "site", values_from = "percent_whole_sp_dry") 
  # arrange(-percent_whole_sp_dry)

table_df %>% 
  gt() %>% 
  data_color(columns = 3,
             colors = scales::col_numeric(palette = gradient_palette, 
                                          domain = c(max(table_df[, 3]), 0))) %>% 
  data_color(columns = 4,
             colors = scales::col_numeric(palette = gradient_palette, 
                                          domain = c(max(table_df[, 4]), 0))) %>% 
  data_color(columns = 5,
             colors = scales::col_numeric(palette = gradient_palette, 
                                          domain = c(max(table_df[, 5]), 0))) %>% 
  data_color(columns = 6,
             colors = scales::col_numeric(palette = gradient_palette, 
                                          domain = c(max(table_df[, 6]), 0))) %>% 
  data_color(columns = 7,
             colors = scales::col_numeric(palette = gradient_palette, 
                                          domain = c(max(table_df[, 7]), 0))) %>% 
  gtsave("sp_site_table-LTE_sites-selection.png", path = here::here("tables"))
