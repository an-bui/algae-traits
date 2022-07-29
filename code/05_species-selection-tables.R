###### 1. overall site table ######

table_df <- biomass %>% 
  filter(year == 2021) %>% 
  select(site, scientific_name, dry_gm2, wm_gm2) %>% 
  filter(site %in% c("aque")) %>% 
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
  arrange(-percent_whole_sp_dry) %>% 
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

###### 2. who is where? ######

# main question: where can I find these species on the transects?
transect_spp <- biomass %>% 
  dplyr::select(-sample_ID) %>% 
  unite("sample_ID", site, year, remove = FALSE) %>% 
  unite("transect_ID", site, year, transect, remove = FALSE) %>% 
  filter(sp_code %in% c(algae_proposal, "LAFA")) %>% 
  filter(site %in% sites_proposal) %>% 
  filter(year == 2021) %>% 
  # shorten CC scientific name
  mutate(scientific_name = case_when(
    sp_code == "CC" ~ "Chondracanthus spp.",
    TRUE ~ scientific_name
  )) %>% 
  mutate(sp_code = dplyr::recode(sp_code, "PH" = "PTCA")) %>% 
  mutate(full = paste(sp_code, " (", scientific_name, ")", sep = "")) %>% 
  select(site, transect, full, dry_gm2) %>% 
  pivot_wider(names_from = full, values_from = dry_gm2)

percov_tbl <- function(site) {
  site_full <- pluck(sites_full, site)
  
  transect_spp %>% 
    filter(site == {{ site }}) %>% 
    select(-site) %>% 
    gt(groupname_col = "transect") %>% 
    tab_header(
      title = site_full,
    ) %>% 
    tab_stubhead(label = "transect") %>% 
    tab_options(
      row_group.as_column = TRUE
    ) %>% 
    # center and middle align transect label
    tab_style(
      style = list(
        cell_text(align = "center", v_align = "middle")
      ),
      locations = cells_row_groups(groups = everything())
    ) 
}

percov_tbl("aque")
percov_tbl("bull")
percov_tbl("napl")
percov_tbl("ivee")
percov_tbl("mohk")

