###### 1. overall site table ######

table_df <- biomass %>%
  filter(between(year, 2011, 2021)) %>%
  select(site, scientific_name, dry_gm2) %>%
  # filter(site %in% c("mohk")) %>%
  # mutate(site = fct_relevel(site, "aque", "napl", "ivee", "mohk", "carp")) %>%
  # filter out MAPY
  filter(scientific_name != "Macrocystis pyrifera",
         site %in% c("bull", "aque", "napl", "ivee", "mohk", "carp")) %>% 
  # create a new column where total biomass across sites for the whole survey is calculated
  # group_by(site) %>%
  mutate(total_dry = sum(dry_gm2)# ,
         # total_wet = sum(wm_gm2)
         ) %>%
  # ungroup() %>%
  # create a new column where total biomass for the species for the whole survey is calculated
  group_by(scientific_name) %>%
  reframe(total_sp_dry = sum(dry_gm2),
         # total_sp_wet = sum(wm_gm2),
         percent_sp_dry = total_sp_dry/total_dry,
         #  # percent_sp_wet = total_sp_wet/total_wet
         ) %>%
  # ungroup() %>%
  unique() %>%
  select(scientific_name, percent_sp_dry) %>%
  mutate(percent_sp_dry = round(percent_sp_dry, 4)) %>%
  mutate(percent_whole_sp_dry = percent_sp_dry*100) %>%
  select(-percent_sp_dry) %>%
  # left_join(., coarse_traits, by = "scientific_name") %>%
  # filter(taxon_phylum == "Rhodophyta")
  # select(scientific_name, sp_code, percent_whole_sp_dry) %>%
  # pivot_wider(names_from = "site", values_from = "percent_whole_sp_dry")
  arrange(-percent_whole_sp_dry)
# 
# table_df %>% 
#   # arrange(-percent_whole_sp_dry) %>% 
#   gt() %>% 
#   data_color(columns = 3,
#              colors = scales::col_numeric(palette = gradient_palette, 
#                                           domain = c(max(table_df[, 3]), 0))) %>% 
#   data_color(columns = 4,
#              colors = scales::col_numeric(palette = gradient_palette, 
#                                           domain = c(max(table_df[, 4]), 0))) %>% 
#   data_color(columns = 5,
#              colors = scales::col_numeric(palette = gradient_palette, 
#                                           domain = c(max(table_df[, 5]), 0))) %>% 
#   data_color(columns = 6,
#              colors = scales::col_numeric(palette = gradient_palette, 
#                                           domain = c(max(table_df[, 6]), 0))) %>% 
#   data_color(columns = 7,
#              colors = scales::col_numeric(palette = gradient_palette, 
#                                           domain = c(max(table_df[, 7]), 0))) %>% 
#   gtsave("sp_site_table-LTE_sites-selection.png", path = here::here("tables"))

###### 2. who is where? ######

# main question: where can I find these species on the transects?
transect_spp <- biomass %>% 
  dplyr::select(-sample_ID) %>% 
  unite("sample_ID", site, year, remove = FALSE) %>% 
  unite("transect_ID", site, year, transect, remove = FALSE) %>% 
  # filter(sp_code %in% c(algae_proposal, "LAFA")) %>% 
  filter(site %in% sites_proposal) %>% 
  filter(year == 2022) %>% 
  # shorten CC scientific name
  mutate(scientific_name = case_when(
    sp_code == "CC" ~ "Chondracanthus spp.",
    sp_code == "DP" ~ "Dictyota spp.", 
    TRUE ~ scientific_name
  )) %>% 
  mutate(sp_code = dplyr::recode(sp_code, "PH" = "PTCA")) %>% 
  mutate(full = paste(sp_code, " (", scientific_name, ")", sep = "")) %>% 
  select(site, transect, full, dry_gm2) %>% 
  mutate(dry_gm2 = round(dry_gm2, 3)) %>% 
  pivot_wider(names_from = full, values_from = dry_gm2) %>% 
  nest(.by = site) %>% 
  arrange(site) %>% 
  # create a new column with full site name
  mutate(site_full = map(site, ~pluck(sites_full, .x))) %>% 
  # make a table using flextable
  mutate(table = map2(data, site_full, ~ 
                        .x %>% 
                        arrange(transect) %>% 
                        flextable() %>% 
                        # making the title of the table the full site name
                        add_header_lines(values = c(.y)) %>% 
                        # aligning everything to be centered
                        align(align = c("center"), part = "header") %>% 
                        align(align = c("center"), part = "body"))) %>% 
  mutate(filename = map(site, ~
                          paste0(.x, "_spp-biomass_", today(), ".docx") %>% 
                          as.vector()))

# saving tables as docx
# map2(transect_spp$table, transect_spp$filename, ~
#          save_as_docx(.x, path = here("tables", "spp-biomass-at-transect", .y)))

transect_rare <- biomass %>% 
  dplyr::select(-sample_ID) %>% 
  unite("sample_ID", site, year, remove = FALSE) %>% 
  unite("transect_ID", site, year, transect, remove = FALSE) %>% 
  filter(!(sp_code %in% c(algae_proposal, "LAFA", "MAPY"))) %>% 
  filter(site %in% sites_proposal) %>% 
  filter(year == 2022) %>% 
  mutate(full = paste(sp_code, " (", scientific_name, ")", sep = "")) %>% 
  select(site, transect, full, dry_gm2) %>% 
  mutate(dry_gm2 = round(dry_gm2, 3)) %>% 
  pivot_wider(names_from = full, values_from = dry_gm2) %>% 
  # nested data frame
  nest(.by = site) %>% 
  arrange(site) %>% 
  # create a new column with full site name
  mutate(site_full = map(site, ~pluck(sites_full, .x))) %>% 
  # select species to keep: only the ones that do occur at the site
  mutate(keepspp = map(data, ~
                         .x %>% 
                         select(!all_of(c("transect"))) %>% 
                         colSums(na.rm = FALSE) %>% 
                         enframe() %>% 
                         filter(value > 0) %>% 
                         pull(name))) %>% 
  # create a data frame to make into a table: only select transect and species to keep
  mutate(table_df = map2(data, keepspp, ~
                           .x %>% 
                           select(transect, all_of(.y)) %>% 
                           arrange(transect))) %>% 
  # make a table using flextable
  mutate(table = map2(table_df, site_full, ~ 
                        .x %>% 
                        flextable() %>% 
                        # making the title of the table the full site name
                        add_header_lines(values = c(.y)) %>% 
                        # aligning everything to be centered
                        align(align = c("center"), part = "header") %>% 
                        align(align = c("center"), part = "body"))) %>% 
  mutate(filename = map(site, ~
                          paste0(.x, "_rare-biomass_", today(), ".docx") %>% 
                          as.vector()))

# saving tables as docx
# map2(transect_rare$table, transect_rare$filename, ~
#          save_as_docx(.x, path = here("tables", "spp-biomass-at-transect", .y)))

# 
# transect_rare %>% 
#   # filter(site == "aque") %>% 
#   flextable() %>% 
#   add_header_row(values = site)
# 
# biomass_tbl <- function(site) {
#   site_full <- pluck(sites_full, site)
#   
#   transect_spp %>% 
#     filter(site == {{ site }}) %>% 
#     select(-site) %>% 
#     gt(groupname_col = "transect") %>% 
#     tab_header(
#       title = site_full,
#     ) %>% 
#     tab_stubhead(label = "transect") %>% 
#     tab_options(
#       row_group.as_column = TRUE
#     ) %>% 
#     # center and middle align transect label
#     tab_style(
#       style = list(
#         cell_text(align = "center", v_align = "middle")
#       ),
#       locations = cells_row_groups(groups = everything())
#     ) 
# }
# 
# biomass_rare_tbl <- function(site) {
#   site_full <- pluck(sites_full, site)
#   
#   transect_rare %>% 
#     filter(site == {{ site }}) %>% 
#     select(-site) %>% 
#     gt(groupname_col = "transect") %>% 
#     tab_header(
#       title = site_full,
#     ) %>% 
#     tab_stubhead(label = "transect") %>% 
#     tab_options(
#       row_group.as_column = TRUE
#     ) %>% 
#     # center and middle align transect label
#     tab_style(
#       style = list(
#         cell_text(align = "center", v_align = "middle")
#       ),
#       locations = cells_row_groups(groups = everything())
#     ) 
# }
# 
# biomass_tbl("bull")
# biomass_tbl("aque")
# biomass_tbl("napl")
# biomass_tbl("ivee")
# biomass_tbl("mohk")
# biomass_tbl("carp")
# 
# biomass_rare_tbl("bull")
# biomass_rare_tbl("aque")
# biomass_rare_tbl("napl")
# biomass_rare_tbl("ivee")
# biomass_rare_tbl("mohk")
# biomass_rare_tbl("carp")

# table_df_percov <- percov %>% 
#   filter(year == 2021) %>% 
#   # select(site, scientific_name, dry_gm2, wm_gm2) %>% 
#   filter(site %in% c("carp")) %>% 
#   # mutate(site = fct_relevel(site, "aque", "napl", "ivee", "mohk", "carp")) %>% 
#   # filter out MAPY
#   filter(scientific_name != "Macrocystis pyrifera") %>% 
#   # create a new column where total biomass for the species for the whole survey is calculated
#   group_by(scientific_name) %>%
#   summarize(sp_percov = sum(percent_cover, na.rm = TRUE)) %>% 
#   mutate(total_percov = sum(sp_percov)) %>% 
#   ungroup() %>% 
#   mutate(sp_prop = round(sp_percov/total_percov, 4)) %>% 
#   mutate(sp_percent = sp_prop*100) %>% 
#   select(-sp_prop) %>% 
#   left_join(., coarse_traits, by = "scientific_name") 

