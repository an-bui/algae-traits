
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 0. source -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# only need to run this once per session
source(here::here("code", "03a_correlation-and-tradeoffs.R"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -------------------- 1. functions and aesthetics ------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# function to extract cumulative proportion, conserved pairwise differences,
# and traits
keep_trait_function <- function(df) {
  df %>% 
    select(contains("trait"), cumu_prop, pairwise_conserved) %>% 
    unnest(cols = c(cumu_prop, pairwise_conserved)) %>% 
    rownames_to_column("combo") %>% 
    mutate(across(contains("trait"), 
                  \(x) fct_relevel(x, trait_factor))) 
}

# function to organize the data frame for the heat map (bottom of upset plot)
keep_traits_heatmap_function <- function(df) {
  df %>% 
    # arrange in descending order of cumulative proportion explained
    arrange(-cumu_prop) %>% 
    # make the order of combinations the order in which they appear
    mutate(combo = fct_inorder((combo))) %>% 
    # get rid of the cumulative proportion column
    select(!cumu_prop) %>% 
    # make the data frame longer for plotting
    pivot_longer(cols = contains("trait")) %>% 
    # rename the colume "value" to be trait
    rename(trait = value) %>% 
    # mark the trait as present
    mutate(present = "TRUE") %>% 
    # get rid of the "name" column, which contains trait1, trait2, etc.
    select(!name) %>% 
    # complete all traits within a combination
    complete(combo, trait) %>% 
    # if the trait is not present in the combination, mark "FALSE"
    replace_na(list(present = "FALSE")) %>% 
    # reorder the traits in factor order
    mutate(trait = fct_relevel(trait, rev(trait_factor)))
}

# aesthetics for the bottom of the upset plot
upset_plot_bottom <- list(
  geom_tile(color = "black"),
  labs(y = "Trait"),
  theme_minimal(),
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        panel.grid = element_blank(),
        text = element_text(size = 14))
)

# aesthetics for the top of the upset plot
upset_plot_top <- list(
  coord_cartesian(ylim = c(0.8, 1)),
  labs(y = "Cumulative proportion"),
  theme_minimal(),
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        panel.grid = element_blank(),
        text = element_text(size = 14))
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------- 2. trait combinations -------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. toy model ----------------------------------------------------------

test_meta <- sample(c(1, 2, 3), replace = TRUE, 20) %>% 
  as_tibble() %>% 
  rename("species" = "value")

test_df <- test_meta %>% 
  mutate(a = case_when(
    species == 1 ~ rnorm(n = 20, mean = 5, sd = 1),
    species == 2 ~ rnorm(n = 20, mean = 10, sd = 1),
    species == 3 ~ rnorm(n = 20, mean = 15, sd = 1)
  )) %>% 
  mutate(b = case_when(
    species == 1 ~ rnorm(n = 20, mean = 6, sd = 1),
    species == 2 ~ rnorm(n = 20, mean = 11, sd = 1),
    species == 3 ~ rnorm(n = 20, mean = 16, sd = 1)
  )) %>% 
  mutate(c = rnorm(n = 20, mean = 8, sd = 1)) %>% 
  mutate(d = rnorm(n = 20, mean = 15, sd = 1)) %>% 
  mutate(e = rnorm(n = 20, mean = 12, sd = 1)) %>% 
  mutate(f = rnorm(n = 20, mean = 7, sd = 1)) %>% 
  mutate(g = rnorm(n = 20, mean = 9, sd = 1)) %>% 
  select(!species)

trait_names <- c("a", "b", "c", "d", "e", "f", "g")

combo <- combn(x = trait_names,
               m = 4) %>% 
  # transpose this to look more like a long format data frame
  t() %>% 
  # turn it into a dataframe
  as_tibble() %>% 
  # rename columns to reflect "traits"
  rename("trait1" = "V1",
         "trait2" = "V2",
         "trait3" = "V3",
         "trait4" = "V4") %>% 
  # nest the data frame
  nest(.by = c(trait1, trait2, trait3, trait4),
       data = everything()) %>% 
  # take out the "data" column (which is meaningless)
  select(!data) %>% 
  # attach the trait data frame to the nested data frame
  # each "cell" contains the trait data frame
  mutate(df = map(
    trait1,
    ~ bind_cols(test_df)
  )) %>% 
  # subset the trait data frame by the traits in the combination
  mutate(subset_df = pmap(
    list(v = df, w = trait1, x = trait2, y = trait3, z = trait4),
    function(v, w, x, y, z) v %>% select(c(w, x, y, z))
  )) %>% 
  # do the PCA
  mutate(pca = map(
    subset_df,
    ~ prcomp(.x, center = TRUE, scale = TRUE)
  )) %>% 
  # extract the cumulative proportion explained by PC1 and PC2
  mutate(cumu_prop = map(
    pca,
    # [3, 2] is cumulative proportion of PC1 and PC2
    ~ summary(.x)$importance[3, 2]
  )) %>% 
  # do the permanova to ask if species are different in trait values
  mutate(permanova = map(
    subset_df,
    ~ adonis2(.x ~ species, data = test_meta)
  )) %>% 
  # do the pairwise comparisons to ask which pairwise differences exist
  mutate(pairwise = map(
    subset_df,
    ~ pairwise.perm.manova(resp = .x, 
                           fact = test_meta$species,
                           p.method = "none")
  )) %>% 
  # ask if the p-values in the pairwise comparisons are greater or less than 0.05
  # if less than 0.05, then subset traits conserve differences between species
  mutate(pairwise_significant = map(
    pairwise,
    ~ .x$p.value < 0.05
  )) %>% 
  # creating a new column: if any p-value > 0.05, then "no" 
  # this makes the following filtering step easier
  mutate(pairwise_conserved = map(
    pairwise_significant,
    ~ case_when(
      # if any p-values < 0.05, then pairwise comparisons are conserved
      "TRUE" %in% .x ~ "yes",
      # if all p-values > 0.05, then pairwise comparisons are not conserved
      TRUE ~ "no"
    )
  ))

# only keep combinations where the pairwise comparisons are conserved
keep <- combo %>% 
  filter(pairwise_conserved == "yes") %>% 
  select(trait1, trait2, trait3, trait4, cumu_prop) %>% 
  unnest(cols = c(cumu_prop))


# ⟞ b. real model ---------------------------------------------------------

# ⟞ ⟞ i. 3 traits ---------------------------------------------------------

# put all the traits into a vector
trait_names_vector <- colnames(pca_mat_log)

# test differences between species with the trait combinations
# first, select 3 traits from the vector (order doesn't matter)
combo_3traits <- combn(x = trait_names_vector,
                       m = 3) %>% 
  # transpose this to look more like a long format data frame
  t() %>% 
  # turn it into a dataframe
  as_tibble() %>% 
  # rename columns to reflect "traits"
  rename("trait1" = "V1",
         "trait2" = "V2",
         "trait3" = "V3") %>% 
  # nest the data frame
  nest(.by = c(trait1, trait2, trait3),
       data = everything()) %>% 
  # take out the "data" column (which is meaningless)
  select(!data) %>% 
  # attach the trait data frame to the nested data frame
  # each "cell" contains the trait data frame
  mutate(df = map(
    trait1,
    ~ bind_cols(pca_mat_log)
  )) %>% 
  # subset the trait data frame by the traits in the combination
  mutate(subset_df = pmap(
    list(v = df, w = trait1, x = trait2, y = trait3),
    function(v, w, x, y) v %>% select(c(w, x, y))
  )) %>% 
  # do the PCA
  mutate(pca = map(
    subset_df,
    ~ prcomp(.x, center = TRUE, scale = TRUE)
  )) %>% 
  # extract the cumulative proportion explained by PC1 and PC2
  mutate(cumu_prop = map(
    pca,
    # [3, 2] is cumulative proportion of PC1 and PC2
    ~ summary(.x)$importance[3, 2]
  )) %>% 
  # do the permanova to ask if species are different in trait values
  mutate(permanova = map(
    subset_df,
    ~ adonis2(.x ~ sp_code, data = ind_traits_filtered)
  )) %>% 
  # do the pairwise comparisons to ask which pairwise differences exist
  mutate(pairwise = map(
    subset_df,
    ~ pairwise.perm.manova(resp = .x, 
                           fact = ind_traits_filtered$sp_code,
                           p.method = "none")
  )) %>% 
  # ask if the p-values in the pairwise comparisons are greater or less than 0.05
  # if less than 0.05, then subset traits conserve differences between species
  mutate(pairwise_significant = map(
    pairwise,
    ~ .x$p.value < 0.05
  )) %>% 
  # creating a new column: if any p-value > 0.05, then "no" 
  # this makes the following filtering step easier
  mutate(pairwise_conserved = map(
    pairwise_significant,
    ~ case_when(
      # if any p-values > 0.05, then pairwise comparisons are NOT conserved
      "FALSE" %in% .x ~ "no",
      # if all p-values > 0.05, then pairwise comparisons are conserved
      TRUE ~ "yes"
    )
  ))

# write_rds(x = combo_3traits,
#           file = here("rds-objects",
#                       paste0("combo_3-traits_", today(), ".rds")))

combo_3traits <- read_rds(file = here("rds-objects",
                                 "combo_3-traits_2024-11-23.rds")) %>% 
  # creating a new column: if any p-value > 0.05, then "no" 
  # this makes the following filtering step easier
  mutate(pairwise_conserved = map(
    pairwise_significant,
    ~ case_when(
      # if any p-values > 0.05, then pairwise comparisons are NOT conserved
      "FALSE" %in% .x ~ "no",
      # if all p-values > 0.05, then pairwise comparisons are conserved
      TRUE ~ "yes"
    )
  ))

# only keep combinations where the pairwise comparisons are conserved
keep_3traits <- combo_3traits %>% 
  keep_trait_function() 

keep_3traits_heatmap <- keep_3traits %>% 
  keep_traits_heatmap_function()

bottom_3traits <- keep_3traits_heatmap %>% 
  ggplot(aes(x = combo,
             y = trait,
             fill = present,
             alpha = pairwise_conserved
  )) +
  scale_fill_manual(values = c("TRUE" = "darkgreen",
                               "FALSE" = "white")) +
  # geom_tile(color = "black")
  scale_alpha_manual(values = c("no" = 0.2,
                                "yes" = 1)) +
  upset_plot_bottom

top_3traits <- keep_3traits %>% 
  # filter(pairwise_conserved == "yes") %>% 
  ggplot(aes(x = reorder(combo, -cumu_prop),
                  y = cumu_prop,
                  alpha = pairwise_conserved)) +
  geom_col(fill = "darkgreen") +
  scale_alpha_manual(values = c("no" = 0.2,
                                "yes" = 1)) +
  upset_plot_top

upset_plot_3traits <- top_3traits / bottom_3traits

upset_plot_3traits

ggsave(here("figures",
            "trait-selection",
            paste0("upset-plot_3traits_all-combinations_", today(), ".jpg")),
       upset_plot_3traits,
       height = 5,
       width = 12,
       units = "cm",
       dpi = 300)

# top combination: SA:V, SA:DW, H:WW
# PCA summary
summary(combo_3traits[[6]][[33]])
# PERMANOVA summary
combo_3traits[[8]][[33]]
# pairwise summary
combo_3traits[[9]][[33]]

# ⟞ ⟞ ii. 4 traits --------------------------------------------------------

combo_4traits <- combn(x = trait_names_vector,
                       m = 4) %>% 
  # transpose this to look more like a long format data frame
  t() %>% 
  # turn it into a dataframe
  as_tibble() %>% 
  # rename columns to reflect "traits"
  rename("trait1" = "V1",
         "trait2" = "V2",
         "trait3" = "V3",
         "trait4" = "V4") %>% 
  # nest the data frame
  nest(.by = c(trait1, trait2, trait3, trait4),
       data = everything()) %>% 
  # take out the "data" column (which is meaningless)
  select(!data) %>% 
  # attach the trait data frame to the nested data frame
  # each "cell" contains the trait data frame
  mutate(df = map(
    trait1,
    ~ bind_cols(pca_mat_log)
  )) %>% 
  # subset the trait data frame by the traits in the combination
  mutate(subset_df = pmap(
    list(v = df, w = trait1, x = trait2, y = trait3, z = trait4),
    function(v, w, x, y, z) v %>% select(c(w, x, y, z))
  )) %>% 
  # do the PCA
  mutate(pca = map(
    subset_df,
    ~ prcomp(.x, center = TRUE, scale = TRUE)
  )) %>% 
  # extract the cumulative proportion explained by PC1 and PC2
  mutate(cumu_prop = map(
    pca,
    # [3, 2] is cumulative proportion of PC1 and PC2
    ~ summary(.x)$importance[3, 2]
  )) %>% 
  # do the permanova to ask if species are different in trait values
  mutate(permanova = map(
    subset_df,
    ~ adonis2(.x ~ sp_code, data = ind_traits_filtered)
  )) %>% 
  # do the pairwise comparisons to ask which pairwise differences exist
  mutate(pairwise = map(
    subset_df,
    ~ pairwise.perm.manova(resp = .x, 
                           fact = ind_traits_filtered$sp_code,
                           p.method = "none")
  )) %>% 
  # ask if the p-values in the pairwise comparisons are greater or less than 0.05
  # if less than 0.05, then subset traits conserve differences between species
  mutate(pairwise_significant = map(
    pairwise,
    ~ .x$p.value < 0.05
  )) %>% 
  # creating a new column: if any p-value > 0.05, then "no" 
  # this makes the following filtering step easier
  mutate(pairwise_conserved = map(
    pairwise_significant,
    ~ case_when(
      # if any p-values > 0.05, then pairwise comparisons are NOT conserved
      "FALSE" %in% .x ~ "no",
      # if all p-values > 0.05, then pairwise comparisons are conserved
      TRUE ~ "yes"
    )
  ))
# write_rds(x = combo_4traits,
#           file = here("rds-objects",
#                       paste0("combo_4-traits_", today(), ".rds")))

combo_4traits <- read_rds(file = here("rds-objects",
                                 "combo_4-traits_2024-11-20.rds")) %>% 
  # creating a new column: if any p-value > 0.05, then "no" 
  # this makes the following filtering step easier
  mutate(pairwise_conserved = map(
    pairwise_significant,
    ~ case_when(
      # if any p-values > 0.05, then pairwise comparisons are NOT conserved
      "FALSE" %in% .x ~ "no",
      # if all p-values > 0.05, then pairwise comparisons are conserved
      TRUE ~ "yes"
    )
  ))

# only keep combinations where the pairwise comparisons are conserved
keep_4traits <- combo_4traits %>% 
  keep_trait_function()

keep_4traits_heatmap <- keep_4traits %>% 
  keep_traits_heatmap_function()

bottom_4traits <- keep_4traits_heatmap %>% 
  ggplot(aes(x = combo,
             y = trait,
             fill = present,
             alpha = pairwise_conserved)) +
  scale_fill_manual(values = c("TRUE" = "cornflowerblue",
                               "FALSE" = "white")) +
  scale_alpha_manual(values = c("no" = 0.2,
                                "yes" = 1)) +
  upset_plot_bottom

top_4traits <- ggplot(data = keep_4traits,
              aes(x = reorder(combo, -cumu_prop),
                  y = cumu_prop,
                  alpha = pairwise_conserved)) +
  geom_col(fill = "cornflowerblue") +
  scale_alpha_manual(values = c("no" = 0.2,
                                "yes" = 1)) +
  upset_plot_top

upset_plot_4traits <- top_4traits / bottom_4traits
upset_plot_4traits

ggsave(here("figures",
            "trait-selection",
            paste0("upset-plot_4traits_all-combinations_", today(), ".jpg")),
       upset_plot_4traits,
       height = 4,
       width = 18,
       units = "cm",
       dpi = 300)

# top combination: SA:V, SA:DW, H, H:WW
# PCA summary
summary(combo_real[[7]][[44]])
# PERMANOVA summary
combo_4traits[[9]][[44]]
# pairwise summary
combo_4traits[[10]][[44]]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ----------------------- 3. reduced model analysis  ----------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. PCA ----------------------------------------------------------------

reduced_mat <- pca_mat_log %>% 
  select(`Surface area`, `Height:wet weight`, `Height`, `Thickness`)

# trait by species PCA
pca_reduced <- rda(reduced_mat, scale = TRUE)

# create a screeplot to visualize axis contributions
screeplot(pca_reduced, bstick = TRUE)

# look at the summary
summary(pca_reduced)

# proportion variance explained for downstream figure making
prop_PC1_reduced <- "66.7%"
prop_PC2_reduced <- "19.7%"

adonis_pairwise_reduced <- pairwise.adonis2(reduced_mat ~ sp_code, 
                                            data = ind_traits_filtered)
rvam_pairwise_reduced <- pairwise.perm.manova(reduced_mat,
                                              fact = ind_traits_filtered$sp_code)
rvam_pairwise_reduced 

# ⟞ b. loadings -----------------------------------------------------------

# get loadings into data frame
loadings_df_reduced <- scores(pca_reduced, 
                              display = 'species', 
                              scaling = 0, 
                              choices = c(1, 2)) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("trait") %>% 
  # arrange the data frame in order of PC1 loadings
  arrange(PC1) %>% 
  # set the factor levels so that all the traits appear in the same order
  mutate(trait = fct_inorder(trait))

pc1_plot_reduced <- ggplot(data = loadings_df_reduced, 
                           aes(x = PC1,
                               y = trait)) +
  loadings_plot_aes +
  scale_x_continuous(limits = c(-1, 1), 
                     breaks = seq(from = -1, to = 1, by = 0.25)) +
  loadings_plot_theme() +
  labs(title = "PC1")

pc2_plot_reduced <- ggplot(data = loadings_df_reduced, 
                           aes(x = PC2,
                               y = trait)) +
  loadings_plot_aes +
  scale_x_continuous(limits = c(-1, 1), 
                     breaks = seq(from = -1, to = 1, by = 0.25)) +
  loadings_plot_theme() +
  labs(title = "PC2")

loadings_plot_reduced <- pc1_plot_reduced / pc2_plot_reduced

ggsave(here::here(
  "figures",
  "ordination",
  paste0("loadings_scale_reduced-model_", today(), ".jpg")),
  loadings_plot_reduced,
  width = 12,
  height = 14,
  units = "cm",
  dpi = 300
)

# ⟞ c. biplots ------------------------------------------------------------

# simple biplot (compare with ggplot output to make sure it's right)
biplot(pca_reduced)

# species points
PCAscores_reduced <- scores(pca_reduced, 
                            display = "sites", 
                            choices = c(1, 2)) %>% 
  as.data.frame() %>% 
  rownames_to_column("specimen_ID") %>% 
  left_join(., metadata_ind, by = "specimen_ID") %>% 
  left_join(., coarse_traits, by = "sp_code") %>% 
  mutate(scientific_name = factor(scientific_name, 
                                  levels = algae_proposal_sciname_factors)) %>% 
  left_join(., fvfm_ind, by = "specimen_ID")

# trait vectors
PCAvect_reduced <- scores(pca_reduced, 
                          display = "species", 
                          choices = c(1, 2)) %>% 
  as.data.frame()

# plot PCA
set.seed(666)
plot_PCA_12_reduced <- ggplot() +
  PCA_aesthetics +
  geom_point(data = PCAscores_reduced, 
             aes(x = PC1, 
                 y = PC2),
             shape = 21,
             color = "darkgrey") +
  geom_segment(data = PCAvect_reduced, 
               aes(x = 0, 
                   y = 0, 
                   xend = PC1, 
                   yend = PC2,
                   color = rownames(PCAvect_reduced)), 
               arrow = arrow(length = unit(0.2, "cm")), 
               linewidth = 1) +
  geom_label_repel(data = PCAvect_reduced, 
                   aes(x = PC1, 
                       y = PC2, 
                       label = rownames(PCAvect_reduced),
                       fill = rownames(PCAvect_reduced)), 
                   size = 8, 
                   alpha = 0.8,
                   color = "black") +
  # scale_x_continuous(limits = c(-2.7, 2.7)) +
  # scale_y_continuous(limits = c(-2.7, 2.7)) +
  PCA_theme() +
  labs(x = paste0("PC1 (", prop_PC1_reduced, ")"),
       y = paste0("PC2 (", prop_PC2_reduced, ")"),
       title = "(a) PC1 and PC2",
       subtitle = "BO, CC, CO, BF, DP, LAFA, PTCA, R, CYOS", 
       color = "Scientific name") 
plot_PCA_12_reduced

ggsave(here::here("figures",
                  "ordination",
                  paste("PCA-log_scale_reduced-model_", today(), ".jpg", sep = "")),
       plot_PCA_12_reduced,
       width = 12, height = 12, units = "cm", dpi = 300)

# ⟞ d. axis contributions -------------------------------------------------

varcoord_reduced <- PCAvect_reduced %>% 
  rownames_to_column("trait") %>% 
  mutate(quality_1 = PC1^2,
         quality_2 = PC2^2) %>% 
  select(trait, quality_1, quality_2) %>% 
  pivot_longer(cols = quality_1:quality_2,
               names_to = "axis",
               values_to = "values") %>% 
  mutate(axis = case_match(
    axis,
    "quality_1" ~ "PC1",
    "quality_2" ~ "PC2"
  )) %>% 
  group_by(axis) %>% 
  mutate(component_total = sum(values)) %>% 
  ungroup() %>% 
  mutate(contrib = (values/component_total)*100)

# The red dashed line on the graph above indicates the expected average contribution. If the contribution of the variables were uniform, the expected value would be 1/length(variables) = 1/10 = 10%
expected_average_reduced <- (1/length(unique(varcoord_reduced$trait)))*100

contrib_aesthetics_reduced <- list(
  geom_col(color = "black"),
  geom_vline(xintercept = expected_average_reduced,
             color = "#96a39e",
             linetype = 2),
  scale_fill_manual(values = trait_color_palette),
  scale_x_continuous(expand = expansion(c(0, 0.05)),
                     position = "top")
)

pc1_contrib_reduced <- varcoord_reduced %>% 
  filter(axis == "PC1") %>% 
  ggplot(aes(y = reorder(trait, contrib),
             x = contrib, 
             fill = trait)) +
  contrib_aesthetics_reduced +
  contrib_theme() +
  labs(title = "(b) Trait % contributions to PC1")

pc1_contrib_reduced

pc2_contrib_reduced <- varcoord_reduced %>% 
  filter(axis == "PC2") %>% 
  ggplot(aes(y = reorder(trait, contrib),
             x = contrib, 
             fill = trait)) +
  contrib_aesthetics_reduced +
  contrib_theme() +
  labs(title = "(c) Trait % contributions to PC2")

pc2_contrib_reduced

# contrib_together_reduced <- pc1_contrib_reduced / pc2_contrib_reduced

contrib_together_reduced <- plot_PCA_12_reduced + (pc1_contrib_reduced / pc2_contrib_reduced)

ggsave(here::here(
  "figures",
  "ordination",
  paste0("contributions_scale_reduced-model_", today(), ".jpg")),
  contrib_together_reduced,
  width = 18,
  height = 10,
  units = "cm",
  dpi = 300
)
