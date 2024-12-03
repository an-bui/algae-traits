# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 0. source -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# only need to run this once per session
source(here::here("code", "03c_trait-combination-visualizations.R"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ----------------------- 1. reduced model analysis  ----------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ⟞ a. 4 traits -----------------------------------------------------------

keep_4trait_combination <- combo_4traits %>% 
  mutate(combo_number = rownames(.)) %>% 
  filter(pairwise_padj_conserved_euc == "yes") %>% 
  unite(trait1:trait4, col = trait_list, sep=", ") %>% 
  mutate(number_of_traits = 4)

keep_3trait_combination <- combo_3traits  %>% 
  mutate(combo_number = rownames(.)) %>% 
  filter(pairwise_padj_conserved_euc == "yes") %>% 
  unite(trait1:trait3, col = trait_list, sep=", ") %>% 
  mutate(number_of_traits = 3)


set.seed(666)

reduced_combinations <- bind_rows(
  keep_4trait_combination,
  keep_3trait_combination
) %>% 
  relocate(number_of_traits, .before = trait_list) %>% 
  relocate(combo_number, .after = number_of_traits) %>% 
  mutate(pca = map(
    subset_df,
    ~ rda(.x, scale = TRUE)
  )) %>%  
  # extract the cumulative proportion explained by PC1 and PC2
  mutate(cumu_prop = map(
    pca,
    # [3, 2] is cumulative proportion of PC1 and PC2
    ~ summary(.x)$cont$importance[3,2]
  )) %>% 
  select(number_of_traits, combo_number, trait_list,
         subset_df, pca, cumu_prop,
         permanova_euc, pairwise_padj_euc) %>% 
  mutate(prop_PC1 = map(
    pca,
    # [2, 1] is proportion explained by PC1
    ~ summary(.x)$cont$importance[2, 1]
  )) %>% 
  mutate(prop_PC2 = map(
    pca,
    # [2, 2] is proportion explained by PC1
    ~ summary(.x)$cont$importance[2, 2] 
  )) %>% 
  mutate(betadisp = map(
    subset_df,
    ~ betadisper(d = vegdist(.x, method = "euclidean"),
                 group = ind_traits_filtered$sp_code)
  )) %>% 
  mutate(betadisp_anova = map(
    betadisp,
    ~ anova(.x)
  )) %>% 
  mutate(betadisp_pairwise = map(
    betadisp,
    ~ TukeyHSD(.x)
  )) %>% 
  mutate(loadings_df = map(
    pca,
    ~ scores(.x, 
             display = 'species', 
             scaling = 0, 
             choices = c(1, 2)) %>% 
      as_tibble(rownames = NA) %>% 
      rownames_to_column("trait") %>% 
      # arrange the data frame in order of PC1 loadings
      arrange(PC1) %>% 
      # set the factor levels so that all the traits appear in the same order
      mutate(trait = fct_inorder(trait))
  )) %>% 
  mutate(pc1_loadings_plot = map(
    loadings_df,
    ~ ggplot(data = .x, 
             aes(x = PC1,
                 y = trait)) +
      loadings_plot_aes +
      scale_x_continuous(limits = c(-1, 1), 
                         breaks = seq(from = -1, to = 1, by = 0.25)) +
      loadings_plot_theme() +
      labs(title = "PC1")
  )) %>% 
  mutate(pc2_loadings_plot = map(
    loadings_df,
    ~ ggplot(data = .x, 
             aes(x = PC2,
                 y = trait)) +
      loadings_plot_aes +
      scale_x_continuous(limits = c(-1, 1), 
                         breaks = seq(from = -1, to = 1, by = 0.25)) +
      loadings_plot_theme() +
      labs(title = "PC2")
  )) %>% 
  mutate(biplot_basic = map(
    pca,
    ~ biplot(.x)
  )) %>% 
  mutate(pca_scores = map(
    pca,
    ~ scores(.x, 
             display = "sites", 
             choices = c(1, 2)) %>% 
      as.data.frame() %>% 
      rownames_to_column("specimen_ID") %>% 
      left_join(., metadata_ind, by = "specimen_ID") %>% 
      left_join(., coarse_traits, by = "sp_code") %>% 
      mutate(scientific_name = factor(scientific_name, 
                                      levels = algae_factors))
  )) %>% 
  mutate(pca_vectors = map(
    pca,
    ~ scores(.x, 
             display = "species", 
             choices = c(1, 2)) %>% 
      as.data.frame()
  )) %>% 
  mutate(biplot_final = pmap(
    list(v = pca_scores, w = pca_vectors, x = prop_PC1, y = prop_PC2),
    function(v, w, x, y) 
      ggplot() +
        PCA_aesthetics +
        vector_colors + 
        geom_point(data = v, 
                   aes(x = PC1, 
                       y = PC2),
                   shape = 21,
                   color = "darkgrey") +
        geom_segment(data = w, 
                     aes(x = 0, 
                         y = 0, 
                         xend = PC1, 
                         yend = PC2,
                         color = rownames(w)), 
                     arrow = arrow(length = unit(0.2, "cm")), 
                     linewidth = 1) +
        geom_label_repel(data = w, 
                         aes(x = PC1, 
                             y = PC2, 
                             label = rownames(w),
                             fill = rownames(w)), 
                         size = 8, 
                         alpha = 0.8,
                         color = "black") +
        # scale_x_continuous(limits = c(-2.7, 2.7)) +
        # scale_y_continuous(limits = c(-2.7, 2.7)) +
        PCA_theme() +
        labs(x = paste0("PC1 (", round(x, 2)*100, "%)"),
             y = paste0("PC2 (", round(y, 2)*100, "%)"),
             title = "(a) PC1 and PC2",
             color = "Scientific name") 
  ))

reduced_combinations[[15]][[1]]
reduced_combinations[[16]][[1]]

reduced_combinations[[20]][[1]]
reduced_combinations[[20]][[2]]
reduced_combinations[[20]][[3]]


# ⟞ ⟞ i. PCA --------------------------------------------------------------

reduced_4traits <- combo_4traits[[6]][[15]]

# trait by species PCA
reduced_4trait_pca <- rda(reduced_4traits, scale = TRUE)

# create a screeplot to visualize axis contributions
screeplot(reduced_4trait_pca, bstick = TRUE)

# look at the summary
summary(reduced_4trait_pca)

# proportion variance explained for downstream figure making
prop_PC1_reduced <- "60%"
prop_PC2_reduced <- "31%"

reduced_permanova <- combo_4traits[[10]][[15]]
reduced_permanova

rvam_pairwise_reduced <- combo_4traits[[14]][[15]]
rvam_pairwise_reduced 

anova(betadisper(d = dist(reduced_mat, method = "euclidean"),
           group = ind_traits_filtered$sp_code))
# different dispersions

TukeyHSD(betadisper(d = dist(reduced_mat, method = "euclidean"),
                    group = ind_traits_filtered$sp_code))
# different dispersions: PTCA-CO, R-BF, PTCA-BF, PTCA-CYOS, PTCA-DP

# ⟞ ⟞ ii. loadings --------------------------------------------------------

# get loadings into data frame
loadings_df_reduced_4traits <- scores(reduced_4trait_pca, 
                                      display = 'species', 
                                      scaling = 0, 
                                      choices = c(1, 2)) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("trait") %>% 
  # arrange the data frame in order of PC1 loadings
  arrange(PC1) %>% 
  # set the factor levels so that all the traits appear in the same order
  mutate(trait = fct_inorder(trait))

pc1_plot_reduced_4traits <- ggplot(data = loadings_df_reduced_4traits, 
                           aes(x = PC1,
                               y = trait)) +
  loadings_plot_aes +
  scale_x_continuous(limits = c(-1, 1), 
                     breaks = seq(from = -1, to = 1, by = 0.25)) +
  loadings_plot_theme() +
  labs(title = "PC1")

pc2_plot_reduced_4traits <- ggplot(data = loadings_df_reduced, 
                           aes(x = PC2,
                               y = trait)) +
  loadings_plot_aes +
  scale_x_continuous(limits = c(-1, 1), 
                     breaks = seq(from = -1, to = 1, by = 0.25)) +
  loadings_plot_theme() +
  labs(title = "PC2")

loadings_plot_reduced_4traits <- pc1_plot_reduced_4traits / pc2_plot_reduced_4traits

# ggsave(here::here(
#   "figures",
#   "ordination",
#   paste0("loadings_scale_reduced_4traits_", today(), ".jpg")),
#   loadings_plot_reduced_4traits,
#   width = 12,
#   height = 14,
#   units = "cm",
#   dpi = 300
# )


# ⟞ ⟞ iii. biplots --------------------------------------------------------

# simple biplot (compare with ggplot output to make sure it's right)
biplot(reduced_4trait_pca)

# species points
PCAscores_reduced_4traits <- scores(reduced_4trait_pca, 
                                    display = "sites", 
                                    choices = c(1, 2)) %>% 
  as.data.frame() %>% 
  rownames_to_column("specimen_ID") %>% 
  left_join(., metadata_ind, by = "specimen_ID") %>% 
  left_join(., coarse_traits, by = "sp_code") %>% 
  mutate(scientific_name = factor(scientific_name, 
                                  levels = algae_factors)) %>% 
  left_join(., fvfm_ind, by = "specimen_ID")

# trait vectors
PCAvect_reduced_4traits <- scores(reduced_4trait_pca, 
                                  display = "species", 
                                  choices = c(1, 2)) %>% 
  as.data.frame()

# plot PCA
set.seed(666)
plot_PCA_12_reduced_4traits <- ggplot() +
  PCA_aesthetics +
  vector_colors + 
  geom_point(data = PCAscores_reduced_4traits, 
             aes(x = PC1, 
                 y = PC2),
             shape = 21,
             color = "darkgrey") +
  geom_segment(data = PCAvect_reduced_4traits, 
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
       color = "Scientific name") 
plot_PCA_12_reduced_4traits

# ggsave(here::here("figures",
#                   "ordination",
#                   paste("PCA-log_scale_4traits_", today(), ".jpg", sep = "")),
#        plot_PCA_12_reduced_4traits,
#        width = 12, height = 12, units = "cm", dpi = 300)


# ⟞ ⟞ iv. axis contributions ----------------------------------------------

varcoord_reduced_4traits <- PCAvect_reduced_4traits %>% 
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
expected_average_reduced_4traits <- (1/length(unique(varcoord_reduced$trait)))*100

contrib_aesthetics_reduced <- list(
  geom_col(color = "black"),
  geom_vline(xintercept = expected_average_reduced_4traits,
             color = "#96a39e",
             linetype = 2),
  scale_fill_manual(values = trait_color_palette),
  scale_x_continuous(expand = expansion(c(0, 0.05)),
                     position = "top")
)

pc1_contrib_reduced_4traits <- varcoord_reduced_4traits %>% 
  filter(axis == "PC1") %>% 
  ggplot(aes(y = reorder(trait, contrib),
             x = contrib, 
             fill = trait)) +
  contrib_aesthetics_reduced +
  contrib_theme() +
  labs(title = "(b) Trait % contributions to PC1")

pc1_contrib_reduced_4traits

pc2_contrib_reduced_4traits <- varcoord_reduced_4traits %>% 
  filter(axis == "PC2") %>% 
  ggplot(aes(y = reorder(trait, contrib),
             x = contrib, 
             fill = trait)) +
  contrib_aesthetics_reduced +
  contrib_theme() +
  labs(title = "(c) Trait % contributions to PC2")

pc2_contrib_reduced

# contrib_together_reduced <- pc1_contrib_reduced / pc2_contrib_reduced

contrib_together_reduced_4traits <- plot_PCA_12_reduced_4traits + (pc1_contrib_reduced_4traits / pc2_contrib_reduced_4traits)

# ggsave(here::here(
#   "figures",
#   "ordination",
#   paste0("contributions_scale_reduced-model_", today(), ".jpg")),
#   contrib_together_reduced,
#   width = 18,
#   height = 10,
#   units = "cm",
#   dpi = 300
# )


# ⟞ b. 3 traits -----------------------------------------------------------

# combination 7: H, T, SA:P
# combination 22: H, D:WW, SA:P

# ⟞ ⟞ i. PCA --------------------------------------------------------------

reduced_3traits_combo7 <- combo_3traits[[5]][[7]]
reduced_3traits_combo22 <- combo_3traits[[5]][[22]]

# trait by species PCA
reduced_3traits_combo7_pca <- rda(reduced_3traits_combo7, scale = TRUE)
reduced_3traits_combo22_pca <- rda(reduced_3traits_combo22, scale = TRUE)

# create a screeplot to visualize axis contributions
screeplot(reduced_3traits_combo7_pca, bstick = TRUE)
screeplot(reduced_3traits_combo22_pca, bstick = TRUE)

# look at the summary
summary(reduced_3traits_combo7_pca)
summary(reduced_3traits_combo22_pca)

# proportion variance explained for downstream figure making
prop_PC1_reduced_3traits_combo7 <- "69%"
prop_PC2_reduced_3traits_combo7 <- "26%"

prop_PC1_reduced_3traits_combo22 <- "78%"
prop_PC2_reduced_3traits_combo22 <- "17%"

reduced_3traits_combo7_permanova <- combo_3traits[[9]][[7]]
reduced_3traits_combo7_permanova

rvam_pairwise_reduced_3traits_combo7 <- combo_3traits[[13]][[7]]
rvam_pairwise_reduced_3traits_combo7

anova(betadisper(d = dist(reduced_3traits_combo7, method = "euclidean"),
                 group = ind_traits_filtered$sp_code))
# different dispersions

TukeyHSD(betadisper(d = dist(reduced_3traits_combo7, method = "euclidean"),
                    group = ind_traits_filtered$sp_code))
# different dispersions: PTCA-CO, PTCA-BO, PTCA-BF, PTCA-CYOS, PTCA-DP

reduced_3traits_combo22_permanova <- combo_3traits[[9]][[22]]
reduced_3traits_combo22_permanova

rvam_pairwise_reduced_3traits_combo22 <- combo_3traits[[13]][[22]]
rvam_pairwise_reduced_3traits_combo22

anova(betadisper(d = dist(reduced_3traits_combo22, method = "euclidean"),
                 group = ind_traits_filtered$sp_code))
# different dispersions

TukeyHSD(betadisper(d = dist(reduced_3traits_combo22, method = "euclidean"),
                    group = ind_traits_filtered$sp_code))
# different dispersions: PTCA-CO, PTCA-BO, PTCA-BF, PTCA-CYOS, PTCA-DP


# ⟞ ⟞ ii. loadings --------------------------------------------------------

# ⟞ ⟞ iii. biplots --------------------------------------------------------

# ⟞ ⟞ iv. axis contributions ----------------------------------------------
