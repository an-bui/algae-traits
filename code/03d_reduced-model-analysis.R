

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ----------------------- 4. reduced model analysis  ----------------------
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

reduced_permanova <- adonis2(reduced_mat ~ sp_code, 
                             data = ind_traits_filtered,
                             method = "euclidean")
reduced_permanova

adonis_pairwise_reduced <- pairwise.adonis2(reduced_mat ~ sp_code, 
                                            data = ind_traits_filtered)
rvam_pairwise_reduced <- pairwise.perm.manova(reduced_mat,
                                              fact = ind_traits_filtered$sp_code,
                                              p.method = "none")
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