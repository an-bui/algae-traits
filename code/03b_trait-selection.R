
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 0. source -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# only need to run this once per session
source(here::here("code", "03a_correlation-and-tradeoffs.R"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------ 1. toy model -----------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ----------------------------- 2. real model -----------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
  mutate(pairwise_padj = map(
    subset_df,
    ~ pairwise.perm.manova(resp = .x, 
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
  # ask if the p-values in the pairwise comparisons are greater or less than 0.05
  # if less than 0.05, then subset traits conserve differences between species
  mutate(pairwise_padj_significant = map(
    pairwise_padj,
    ~ .x$p.value < 0.05
  )) %>% 
  # creating a new column: if any p-value > 0.05, then "no" 
  # this makes the following filtering step easier
  mutate(pairwise_padj_conserved = map(
    pairwise_padj_significant,
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


# ⟞ ii. 3 traits ----------------------------------------------------------

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
  mutate(pairwise_padj = map(
    subset_df,
    ~ pairwise.perm.manova(resp = .x, 
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
  # ask if the p-values in the pairwise comparisons are greater or less than 0.05
  # if less than 0.05, then subset traits conserve differences between species
  mutate(pairwise_padj_significant = map(
    pairwise_padj,
    ~ .x$p.value < 0.05
  )) %>% 
  # creating a new column: if any p-value > 0.05, then "no" 
  # this makes the following filtering step easier
  mutate(pairwise_padj_conserved = map(
    pairwise_padj_significant,
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




