
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 0. source -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# only need to run this once per session
source(here::here("code", "03a_correlation-and-tradeoffs.R"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -------------------------- 1. trait selection  --------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# put all the traits into a vector
trait_names_vector <- colnames(pca_mat_log)

# ⟞ i. 4 traits -----------------------------------------------------------

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
    function(v, w, x, y, z) v %>% select(all_of(c(w, x, y, z)))
  )) %>% 
  # do the PCA
  mutate(pca = map(
    subset_df,
    ~ rda(.x, scale = TRUE)
  )) %>% 
  # extract the cumulative proportion explained by PC1 and PC2
  mutate(cumu_prop = map(
    pca,
    # [3, 2] is cumulative proportion of PC1 and PC2
    ~ summary(.x)$cont$importance[3, 2]
  )) %>% 
  mutate(dist_obj = map(
    subset_df,
    ~ vegdist(.x, method = "euclidean")
  )) %>% 
  # do the permanova to ask if species are different in trait values
  # mutate(permanova = map(
  #   subset_df,
  #   ~ adonis2(.x ~ sp_code, 
  #             data = ind_traits_filtered)
  # )) %>% 
  # repeat permanova with euclidean distances
  mutate(permanova_euc = map(
    dist_obj,
    ~ adonis2(.x ~ sp_code, 
              data = ind_traits_filtered)
  )) %>% 
  # do the pairwise comparisons to ask which pairwise differences exist
  # mutate(pairwise = map(
  #   subset_df,
  #   ~ pairwise.perm.manova(resp = .x, 
  #                          fact = ind_traits_filtered$sp_code,
  #                          p.method = "none")
  # )) %>% 
  # mutate(pairwise_padj = map(
  #   subset_df,
  #   ~ pairwise.perm.manova(resp = .x, 
  #                          fact = ind_traits_filtered$sp_code,
  #                          p.method = "BH")
  # )) %>% 
  mutate(pairwise_euc = map(
    dist_obj,
    ~ pairwise.perm.manova(resp = .x, 
                           fact = ind_traits_filtered$sp_code,
                           p.method = "none")
  )) %>% 
  mutate(pairwise_padj_euc = map(
    dist_obj,
    ~ pairwise.perm.manova(resp = .x, 
                           fact = ind_traits_filtered$sp_code,
                           p.method = "BH")
  )) %>% 
  # ask if the p-values in the pairwise comparisons are greater or less than 0.05
  # if less than 0.05, then subset traits conserve differences between species
  # mutate(pairwise_significant = map(
  #   pairwise,
  #   ~ .x$p.value < 0.05
  # )) %>% 
  # creating a new column: if any p-value > 0.05, then "no" 
  # this makes the following filtering step easier
  # mutate(pairwise_conserved = map(
  #   pairwise_significant,
  #   ~ case_when(
  #     # if any p-values > 0.05, then pairwise comparisons are NOT conserved
  #     "FALSE" %in% .x ~ "no",
  #     # if all p-values > 0.05, then pairwise comparisons are conserved
  #     TRUE ~ "yes"
  #   )
  # )) %>% 
  # # repeat for p-value adjusted (not using euclidean distances)
  # mutate(pairwise_padj_significant = map(
  #   pairwise_padj,
  #   ~ .x$p.value < 0.05
  # )) %>% 
  # mutate(pairwise_padj_conserved = map(
  #   pairwise_padj_significant,
  #   ~ case_when(
  #     "FALSE" %in% .x ~ "no",
  #     TRUE ~ "yes"
  #   )
  # )) %>% 
  mutate(pairwise_significant_euc = map(
    pairwise_euc,
    ~ .x$p.value < 0.05
  )) %>% 
  mutate(match_pairs = map(
    pairwise_significant_euc,
    ~ .x == full_pairwise_matrix
  )) %>% 
  # creating a new column: if any p-value > 0.05, then "no" 
  # this makes the following filtering step easier
  mutate(pairwise_conserved_euc = map(
    match_pairs,
    ~ case_when(
      # if any p-values > 0.05, then pairwise comparisons are NOT conserved
      "FALSE" %in% .x ~ "no",
      # if all p-values > 0.05, then pairwise comparisons are conserved
      TRUE ~ "yes"
    )
  )) %>% 
  mutate(pairwise_padj_significant_euc = map(
    pairwise_padj_euc,
    ~ .x$p.value < 0.05
  )) %>% 
  mutate(pairwise_padj_conserved_euc = map(
    pairwise_padj_significant_euc,
    ~ case_when(
      "FALSE" %in% .x ~ "no",
      TRUE ~ "yes"
    )
  ))

# write_rds(x = combo_4traits,
#           file = here("rds-objects",
#                       paste0("combo_4-traits_", today(), ".rds")))


# ⟞ ii. 3 traits ----------------------------------------------------------

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
    function(v, w, x, y) v %>% select(all_of(c(w, x, y)))
  )) %>% 
  # do the PCA
  mutate(pca = map(
    subset_df,
    ~ rda(.x, scale = TRUE)
  )) %>% 
  # extract the cumulative proportion explained by PC1 and PC2
  mutate(cumu_prop = map(
    pca,
    # [3, 2] is cumulative proportion of PC1 and PC2
    ~ summary(.x)$cont$importance[3, 2]
  )) %>% 
  mutate(dist_obj = map(
    subset_df,
    ~ vegdist(.x, method = "euclidean")
  )) %>% 
  # do the permanova to ask if species are different in trait values
  # mutate(permanova = map(
  #   subset_df,
  #   ~ adonis2(.x ~ sp_code, 
  #             data = ind_traits_filtered)
  # )) %>% 
  # repeat permanova with euclidean distances
  mutate(permanova_euc = map(
    dist_obj,
    ~ adonis2(.x ~ sp_code, 
              data = ind_traits_filtered)
  )) %>% 
  # do the pairwise comparisons to ask which pairwise differences exist
  # mutate(pairwise = map(
  #   subset_df,
  #   ~ pairwise.perm.manova(resp = .x, 
  #                          fact = ind_traits_filtered$sp_code,
  #                          p.method = "none")
  # )) %>% 
  # mutate(pairwise_padj = map(
  #   subset_df,
  #   ~ pairwise.perm.manova(resp = .x, 
  #                          fact = ind_traits_filtered$sp_code,
  #                          p.method = "BH")
  # )) %>% 
  mutate(pairwise_euc = map(
    dist_obj,
    ~ pairwise.perm.manova(resp = .x, 
                           fact = ind_traits_filtered$sp_code,
                           p.method = "none")
  )) %>% 
  mutate(pairwise_padj_euc = map(
    dist_obj,
    ~ pairwise.perm.manova(resp = .x, 
                           fact = ind_traits_filtered$sp_code,
                           p.method = "BH")
  )) %>% 
  # ask if the p-values in the pairwise comparisons are greater or less than 0.05
  # if less than 0.05, then subset traits conserve differences between species
  # mutate(pairwise_significant = map(
  #   pairwise,
  #   ~ .x$p.value < 0.05
  # )) %>% 
  # creating a new column: if any p-value > 0.05, then "no" 
  # this makes the following filtering step easier
  # mutate(pairwise_conserved = map(
  #   pairwise_significant,
  #   ~ case_when(
  #     # if any p-values > 0.05, then pairwise comparisons are NOT conserved
  #     "FALSE" %in% .x ~ "no",
  #     # if all p-values > 0.05, then pairwise comparisons are conserved
  #     TRUE ~ "yes"
  #   )
  # )) %>% 
  # # repeat for p-value adjusted (not using euclidean distances)
  # mutate(pairwise_padj_significant = map(
  #   pairwise_padj,
  #   ~ .x$p.value < 0.05
  # )) %>% 
  # mutate(pairwise_padj_conserved = map(
  #   pairwise_padj_significant,
  #   ~ case_when(
  #     "FALSE" %in% .x ~ "no",
  #     TRUE ~ "yes"
  #   )
  # )) %>% 
  mutate(pairwise_significant_euc = map(
    pairwise_euc,
    ~ .x$p.value < 0.05
  )) %>% 
  mutate(match_pairs = map(
    pairwise_significant_euc,
    ~ .x == full_pairwise_matrix
  )) %>% 
  # creating a new column: if any p-value > 0.05, then "no" 
  # this makes the following filtering step easier
  mutate(pairwise_conserved_euc = map(
    match_pairs,
    ~ case_when(
      # if any p-values > 0.05, then pairwise comparisons are NOT conserved
      "FALSE" %in% .x ~ "no",
      # if all p-values > 0.05, then pairwise comparisons are conserved
      TRUE ~ "yes"
    )
  )) %>% 
  mutate(pairwise_padj_significant_euc = map(
    pairwise_padj_euc,
    ~ .x$p.value < 0.05
  )) %>% 
  mutate(pairwise_padj_conserved_euc = map(
    pairwise_padj_significant_euc,
    ~ case_when(
      "FALSE" %in% .x ~ "no",
      TRUE ~ "yes"
    )
  ))

# write_rds(x = combo_3traits,
#           file = here("rds-objects",
#                       paste0("combo_3-traits_", today(), ".rds")))


# ⟞ iii. 2 traits ---------------------------------------------------------

# test differences between species with the trait combinations
# first, select 2 traits from the vector (order doesn't matter)
combo_2traits <- combn(x = trait_names_vector,
                       m = 2) %>% 
  # transpose this to look more like a long format data frame
  t() %>% 
  # turn it into a dataframe
  as_tibble() %>% 
  # rename columns to reflect "traits"
  rename("trait1" = "V1",
         "trait2" = "V2") %>% 
  # nest the data frame
  nest(.by = c(trait1, trait2),
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
    list(v = df, w = trait1, x = trait2),
    function(v, w, x) v %>% select(all_of(c(w, x)))
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
    ~ adonis2(.x ~ sp_code, 
              data = ind_traits_filtered)
  )) %>% 
  # repeat permanova with euclidean distances
  mutate(permanova_euc = map(
    subset_df,
    ~ adonis2(.x ~ sp_code, 
              data = ind_traits_filtered, 
              method = "euclidean")
  )) %>% 
  # do the pairwise comparisons to ask which pairwise differences exist
  mutate(pairwise = map(
    subset_df,
    ~ pairwise.perm.manova(resp = .x, 
                           fact = ind_traits_filtered$sp_code,
                           p.method = "none")
  )) %>% 
  mutate(pairwise_padj = map(
    subset_df,
    ~ pairwise.perm.manova(resp = .x, 
                           fact = ind_traits_filtered$sp_code,
                           p.method = "BH")
  )) %>% 
  mutate(pairwise_euc = map(
    subset_df,
    ~ pairwise.perm.manova(resp = dist(.x, "euclidean"), 
                           fact = ind_traits_filtered$sp_code,
                           p.method = "none")
  )) %>% 
  mutate(pairwise_padj_euc = map(
    subset_df,
    ~ pairwise.perm.manova(resp = dist(.x, "euclidean"), 
                           fact = ind_traits_filtered$sp_code,
                           p.method = "BH")
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
  )) %>% 
  # repeat for p-value adjusted (not using euclidean distances)
  mutate(pairwise_padj_significant = map(
    pairwise_padj,
    ~ .x$p.value < 0.05
  )) %>% 
  mutate(pairwise_padj_conserved = map(
    pairwise_padj_significant,
    ~ case_when(
      "FALSE" %in% .x ~ "no",
      TRUE ~ "yes"
    )
  )) %>% 
  mutate(pairwise_significant_euc = map(
    pairwise_euc,
    ~ .x$p.value < 0.05
  )) %>% 
  # creating a new column: if any p-value > 0.05, then "no" 
  # this makes the following filtering step easier
  mutate(pairwise_conserved_euc = map(
    pairwise_significant_euc,
    ~ case_when(
      # if any p-values > 0.05, then pairwise comparisons are NOT conserved
      "FALSE" %in% .x ~ "no",
      # if all p-values > 0.05, then pairwise comparisons are conserved
      TRUE ~ "yes"
    )
  )) %>% 
  mutate(pairwise_padj_significant_euc = map(
    pairwise_padj_euc,
    ~ .x$p.value < 0.05
  )) %>% 
  mutate(pairwise_padj_conserved_euc = map(
    pairwise_padj_significant_euc,
    ~ case_when(
      "FALSE" %in% .x ~ "no",
      TRUE ~ "yes"
    )
  ))

# write_rds(x = combo_2traits,
#           file = here("rds-objects",
#                       paste0("combo_2-traits_", today(), ".rds")))




# ⟞ iv. 5 traits ----------------------------------------------------------

# test differences between species with the trait combinations
# first, select 3 traits from the vector (order doesn't matter)
combo_5traits <- combn(x = trait_names_vector,
                       m = 5) %>% 
  # transpose this to look more like a long format data frame
  t() %>% 
  # turn it into a dataframe
  as_tibble() %>% 
  # rename columns to reflect "traits"
  rename("trait1" = "V1",
         "trait2" = "V2",
         "trait3" = "V3",
         "trait4" = "V4",
         "trait5" = "V5") %>% 
  # nest the data frame
  nest(.by = c(trait1, trait2, trait3, trait4, trait5),
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
    list(u = df, v = trait1, w = trait2, x = trait3, y = trait4, z = trait5),
    function(u, v, w, x, y, z) u %>% select(all_of(c(v, w, x, y, z)))
  )) %>% 
  # do the PCA
  mutate(pca = map(
    subset_df,
    ~ rda(.x, scale = TRUE)
  )) %>% 
  # extract the cumulative proportion explained by PC1 and PC2
  mutate(cumu_prop = map(
    pca,
    # [3, 2] is cumulative proportion of PC1 and PC2
    ~ summary(.x)$cont$importance[3, 2]
  )) %>% 
  # create a distance object
  mutate(dist_obj = map(
    subset_df,
    ~ vegdist(.x, method = "euclidean")
  )) %>% 
  # do the permanova to ask if species are different in trait values
  mutate(permanova_euc = map(
    dist_obj,
    ~ adonis2(.x ~ sp_code, 
              data = ind_traits_filtered)
  )) %>% 
  # do the pairwise comparisons to ask which pairwise differences exist
  mutate(pairwise_euc = map(
    dist_obj,
    ~ pairwise.perm.manova(resp = .x, 
                           fact = ind_traits_filtered$sp_code,
                           p.method = "none")
  )) %>% 
  mutate(pairwise_padj_euc = map(
    dist_obj,
    ~ pairwise.perm.manova(resp = .x, 
                           fact = ind_traits_filtered$sp_code,
                           p.method = "BH")
  )) %>% 
  # ask if the p-values in the pairwise comparisons are greater or less than 0.05
  # if less than 0.05, then subset traits conserve differences between species
  mutate(pairwise_significant_euc = map(
    pairwise_euc,
    ~ .x$p.value < 0.05
  )) %>% 
  # creating a new column: if any p-value > 0.05, then "no" 
  # this makes the following filtering step easier
  mutate(pairwise_conserved_euc = map(
    pairwise_significant_euc,
    ~ case_when(
      # if any p-values > 0.05, then pairwise comparisons are NOT conserved
      "FALSE" %in% .x ~ "no",
      # if all p-values > 0.05, then pairwise comparisons are conserved
      TRUE ~ "yes"
    )
  )) %>% 
  # repeating for adjusted p-values
  mutate(pairwise_padj_significant_euc = map(
    pairwise_padj_euc,
    ~ .x$p.value < 0.05
  )) %>% 
  mutate(pairwise_padj_conserved_euc = map(
    pairwise_padj_significant_euc,
    ~ case_when(
      "FALSE" %in% .x ~ "no",
      TRUE ~ "yes"
    )
  ))

# no combination of 5 traits yields all significant pairwise differences

# write_rds(x = combo_5traits,
#           file = here("rds-objects",
#                       paste0("combo_5-traits_", today(), ".rds")))
