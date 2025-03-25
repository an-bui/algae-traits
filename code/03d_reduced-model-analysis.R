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
  # count TRUE and FALSE occurrences
  mutate(count = map(
    pairwise_significant_euc,
    ~ .x %>% unlist() %>% table() 
  )) %>% 
  # number of insignificant pairwise comparisons
  mutate(count_false = pmap(
    list(x = pairwise_padj_conserved_euc, y = count),
    function(x, y) case_when(
      x == "no" ~ y[[1]],
      x == "yes" ~ 0
    )
  )) %>% 
  filter(count_false %in% c(0)) %>% 
  # filter(pairwise_padj_conserved_euc == "yes") %>% 
  unite(trait1:trait4, col = trait_list, sep=", ") %>% 
  mutate(number_of_traits = 4)

keep_3trait_combination <- combo_3traits %>% 
  mutate(combo_number = rownames(.)) %>% 
  # count TRUE and FALSE occurrences
  mutate(count = map(
    pairwise_padj_significant_euc,
    ~ .x %>% unlist() %>% table() 
  )) %>% 
  # number of insignificant pairwise comparisons
  mutate(count_false = pmap(
    list(x = pairwise_padj_conserved_euc, y = count),
    function(x, y) case_when(
      x == "no" ~ y[[1]],
      x == "yes" ~ 0
    )
  )) %>% 
  filter(count_false %in% c(0)) %>% 
  filter(!(combo_number == 3)) %>% 
  unite(trait1:trait3, col = trait_list, sep=", ") %>% 
  mutate(number_of_traits = 3)


set.seed(666)

reduced_combinations <- bind_rows(
  keep_4trait_combination,
  keep_3trait_combination
) %>% 
  relocate(number_of_traits, .before = trait_list) %>% 
  relocate(combo_number, .after = number_of_traits) %>% 
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
                                      levels = algae_factors)) %>% 
      mutate(sp_label = case_when(
        scientific_name == "Dictyota binghamiae; Dictyota flabellata; Dictyota coriacea" ~ "Dictyota spp.",
        scientific_name == "Chondracanthus corymbiferus; Chondracanthus exasperatus" ~ "Chondracanthus spp.",
        TRUE ~ scientific_name
      )) %>% 
      mutate(sp_label = fct_relevel(sp_label, algae_splabel_factors))
  )) %>% 
  mutate(pca_vectors = map(
    pca,
    ~ scores(.x, 
             display = "species", 
             choices = c(1, 2)) %>% 
      as.data.frame()
  )) %>% 
  mutate(vector_biplot_final = pmap(
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
  )) %>% 
  mutate(species_biplot_final = pmap(
    list(v = pca_scores, x = prop_PC1, y = prop_PC2),
    function(v, x, y)
    ggplot(data = v, 
           aes(x = PC1, 
               y = PC2, 
               color = sp_label, 
    ) ) +
      PCA_aesthetics +
      species_colors + 
      geom_point(
        shape = 21,
        alpha = 0.3,
        size = 1
      )   +
      stat_ellipse(aes(color = sp_label),
                   level = 0.5) +
      # scale_x_continuous(limits = c(-1.3, 1.3)) +
      # scale_y_continuous(limits = c(-1.3, 1.3)) +
      PCA_theme() +
      theme(legend.position = "right",
            legend.text = element_text(size = 16),
            legend.background = element_blank(),
            legend.spacing.y = unit(0.01, "cm"),
            legend.key.spacing.y = unit(0.01, "cm"),
            legend.key.height = unit(0.25, "cm")) +
      labs(color = "Scientific name",
           title = "(b) Species",
           x = paste0("PC1 (", round(x, 2)*100, "%)"),
           y = paste0("PC2 (", round(y, 2)*100, "%)"),
      )
  )) %>% 
  mutate(plots_together = pmap(
    list(x = vector_biplot_final, y = species_biplot_final),
    function(x, y) x + y
  ))

# plots together
# 3 traits combo 7: row 69
# 3 traits combo 22: row 73
# 4 traits combo 15: row 8
combo_4traits_combo2_biplots <- reduced_combinations[[22]][[1]]
combo_4traits_combo5_biplots <- reduced_combinations[[22]][[2]]
combo_3traits_combo9_biplots <- reduced_combinations[[22]][[3]]

combo_biplots_together <- list(
  combo_4traits_combo2_biplots,
  combo_4traits_combo5_biplots,
  combo_3traits_combo9_biplots
)

combo_biplot_file_names <- list(
  "biplot-4traits-combo2",
  "biplot-4traits-combo5",
  "biplot-3traits-combo9"
)

for(i in 1:length(combo_biplots_together)) {
  ggsave(here("figures",
              "trait-selection", 
              paste0("reduced-model_", 
                     combo_biplot_file_names[[i]], 
                     "_", 
                     today(), 
                     ".jpg")),
         plot = combo_biplots_together[[i]],
         width = 24,
         height = 12,
         units = "cm",
         dpi = 300)
}

dispersions <- reduced_combinations %>% 
  select(number_of_traits, combo_number, betadisp_anova) %>% 
  mutate(p_val = map(
    betadisp_anova,
    ~ case_when(
      .x$`Pr(>F)`[[1]] > 0.05 ~ paste0("no diff (p = ", round(.x$`Pr(>F)`[[1]], digits = 2), ")"),
      .x$`Pr(>F)`[[1]] < 0.05 ~ paste0("diff dispersions (p = ", round(.x$`Pr(>F)`[[1]], digits = 2), ")")
    )
   # ~ .x$`Pr(>F)`[[1]]
  ))

