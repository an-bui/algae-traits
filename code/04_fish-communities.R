# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 0. source -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# only need to run this once per session
source(here::here("code", "01a_trait-cleaning.R"))



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

fd_feve <- fundiversity::fd_feve(
  traits = group_df,
  sp_com = algae_matrix
)

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
           y = fish_sppnum,
           color = site)) +
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