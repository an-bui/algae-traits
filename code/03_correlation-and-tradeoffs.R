
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 0. source -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# only need to run this once per session
source(here::here("code", "01a_trait-cleaning.R"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------ 1. wrangling -----------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pca_mat <- ind_traits %>% 
  filter(sp_code %in% algae_proposal) %>% 
  # weird Nienburgia?
  filter(!(specimen_ID %in% c("20210721-BULL-023", "20210719-IVEE-009"))) %>% 
  # filter(sp_code != "EGME") %>% 
  mutate(growth_form_num = case_when(
    growth_form == "leathery_macrophyte" ~ 1,
    growth_form == "corticated_macrophytes" ~ 2,
    growth_form == "corticated_foliose" ~ 3
  ), 
  pigment_type_num = case_when(
    pigment_type == "red" ~ 1,
    pigment_type == "brown" ~ 2
  ), 
  longevity_num = case_when(
    longevity == "annual" ~ 1,
    longevity == "perennial" ~ 3,
    longevity == "annual or perennial" ~ 2
  )) %>% 
  select(specimen_ID, 
         maximum_height, mass_to_height, sav_mean, thickness_mm_mean, 
         tdmc_mean, sta_mean, sav_mean, sap_mean, fvfm_mean, 
         aspect_ratio_mean) %>% 
  column_to_rownames("specimen_ID") %>% 
  drop_na() %>% 
  rename(`Maximum height` = maximum_height,
         `Mass:height` = mass_to_height,
         `SA:V` = sav_mean,
         `Thickness` = thickness_mm_mean,
         `TDMC` = tdmc_mean,
         `STA` = sta_mean,
         `SA:V` = sav_mean,
         `SA:P` = sap_mean,
         `Fv/Fm` = fvfm_mean,
         `Aspect ratio` = aspect_ratio_mean)

pca_mat_scale <- scale(pca_mat)

pca_mat_log <- pca_mat %>% 
  # only log transforming traits that were not normally distributed
  mutate(across(c(`Maximum height`, `Mass:height`, `SA:V`,
                  `Thickness`, `STA`, `SA:P`, `Aspect ratio`), 
                log))

pca_mat_with_cn <- ind_traits %>% 
  filter(sp_code %in% algae_proposal) %>% 
  # weird Nienburgia?
  filter(!(specimen_ID %in% c("20210721-BULL-023", "20210719-IVEE-009"))) %>% 
  mutate(growth_form_num = case_when(
    growth_form == "leathery_macrophyte" ~ 1,
    growth_form == "corticated_macrophytes" ~ 2,
    growth_form == "corticated_foliose" ~ 3
  ), 
  pigment_type_num = case_when(
    pigment_type == "red" ~ 1,
    pigment_type == "brown" ~ 2
  ), 
  longevity_num = case_when(
    longevity == "annual" ~ 1,
    longevity == "perennial" ~ 3,
    longevity == "annual or perennial" ~ 2
  )) %>% 
  select(specimen_ID, 
         maximum_height, total_volume, sav_mean, thickness_mm_mean, tdmc_mean, 
         sta_mean, sav_mean, sap_mean, fvfm_mean, thickness_mm_mean, 
         weight_percent_n_mean, weight_percent_c_mean) %>% 
  column_to_rownames("specimen_ID") %>% 
  drop_na()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ----------------------------- 2. correlations ---------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mat_for_corr <- pca_mat_log

# create a correlation matrix based on pearson correlation
M <- mat_for_corr %>% 
  cor(method = "pearson")

p <- cor.mtest(mat_for_corr,
               conf.level = 0.95,
               method = "pearson")

corr_traits <- corrplot(corr = M, 
                        p.mat = p$p,
                        diag = FALSE,
                        method = "circle", 
                        # addgrid.col = NA,
                        type = "lower", 
                        insig = "blank", 
                        addCoef.col = "black", 
                        col= colorRampPalette(c("#F21A00", "#FFFFFF", "#3B9AB2"))(200), 
                        sig.level = 0.05)

jpeg(here::here("figures", 
                "correlation",
                paste0("corrplot_log_", today(), ".jpeg")),
     width = 14, height = 14, units = "cm", res = 200)
corrplot(corr = M, 
         p.mat = p$p,
         diag = FALSE,
         method = "circle", 
         # addgrid.col = NA,
         type = "lower", 
         insig = "blank", 
         addCoef.col = "black", 
         col= colorRampPalette(c("#F21A00", "#FFFFFF", "#3B9AB2"))(200), 
         sig.level = 0.05)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------ 3. standardized major axis regression ----------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. identify significant correlations ----------------------------------

# p-values from correlation in matrix form
p_df <- p$p

# thank you rguha for this solution: https://stackoverflow.com/questions/24554642/what-is-the-best-way-to-remove-items-above-a-diagonal-in-a-data-frame-in-r
# get rid of the upper triangle
p_df[upper.tri(p_df)] <- NA

pairwise_comparisons <- p_df %>% 
  # turn into data frame
  as.data.frame() %>% 
  # make the rownames into "trait1" or the first trait in the pairwise comparison
  rownames_to_column("trait1") %>% 
  # when p < 0.05, keep the p, and when it's not then replace with NA
  mutate(across(where(is.numeric), 
                ~ case_when(. < 0.05 ~ ., TRUE ~ NA))) %>% 
  # make the data frame longer and turn the column names into "trait2" or the second trait in the comparison
  pivot_longer(cols = `Maximum height`:`Aspect ratio`,
               names_to = "trait2",
               values_to = "corr") %>%
  # drop any NAs in the p-values (the upper triangle of the matrix)
  drop_na(corr) %>% 
  # only include p-values that are more than 0 (the matrix has the same-same pair as 0)
  filter(corr > 0) %>% 
  mutate(trait1_column = case_match(
    trait1,
    "Mass:height" ~ "mass_to_height",
    "SA:V" ~ "sav_mean",
    "Thickness" ~ "thickness_mm_mean",
    "TDMC" ~ "tdmc_mean",
    "STA" ~ "sta_mean",
    "SA:P" ~ "sap_mean",
    "Fv/Fm" ~ "fvfm_mean",
    "Aspect ratio" ~ "aspect_ratio_mean"
  )) %>% 
  mutate(trait2_column = case_match(
    trait2,
    "Maximum height" ~ "maximum_height",
    "Mass:height" ~ "mass_to_height",
    "SA:V" ~ "sav_mean",
    "Thickness" ~ "thickness_mm_mean",
    "TDMC" ~ "tdmc_mean",
    "STA" ~ "sta_mean",
    "SA:P" ~ "sap_mean",
    "Fv/Fm" ~ "fvfm_mean"
  )) %>% 
  mutate(formula = paste(trait1_column, "~", trait2_column, sep = " "))

log_ind_traits <- ind_traits %>% 
  mutate(across(c(maximum_height, mass_to_height, sav_mean,
                  thickness_mm_mean, sta_mean, sap_mean, aspect_ratio_mean), 
                log))

scaled_ind_traits <- ind_traits_filtered %>% 
  mutate(across(where(is.numeric), scale))

pairwise_sma <- function(model_formula, trait1, trait2, data) {
  df <- if(data == "normal") {
    ind_traits_filtered
  } else if(data == "scaled") {
    scaled_ind_traits
  } else if(data == "log") {
    log_ind_traits
  } 
  
  # lmodel2
  lmodel2_obj <- lmodel2(model_formula,
                         data = df)
  # sma
  sma_obj <- sma(model_formula,
                 data = df)
  
  # plot
  plot <- ggplot(data = df,
                 aes(x = {{ trait1 }},
                     y = {{ trait2 }})) +
    geom_point(alpha = 0.1) +
    stat_ma_line(method = "SMA") +
    # labs(title = model_formula) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  # return all
  return(list(lmodel2_obj, sma_obj, plot))
  
  
}

pair_sav_height <- pairwise_sma(
  model_formula = "sav_mean ~ maximum_height", 
  trait1 = sav_mean,
  trait2 = maximum_height,
  data = "log"
)

pair_sav_mh <- pairwise_sma(
  model_formula = "sav_mean ~ mass_to_height", 
  trait1 = sav_mean,
  trait2 = mass_to_height,
  data = "log"
)

# fix this later with new mass data entry, though not sure why things worked with 
# other models...
pair_thickness_mh <- pairwise_sma(
  model_formula = "thickness_mm_mean ~ mass_to_height",
  trait1 = thickness_mm_mean,
  trait2 = mass_to_height,
  data = "log"
)

pair_thickness_sav <- pairwise_sma(
  model_formula = "thickness_mm_mean ~ sav_mean", 
  trait1 = thickness_mm_mean,
  trait2 = sav_mean,
  data = "log"
)

pair_tdmc_height <- pairwise_sma(
  model_formula = "tdmc_mean ~ maximum_height", 
  trait1 = tdmc_mean,
  trait2 = maximum_height,
  data = "log"
)

pair_tdmc_mh <- pairwise_sma(
  model_formula = "tdmc_mean ~ mass_to_height", 
  trait1 = tdmc_mean,
  trait2 = mass_to_height,
  data = "log"
)

pair_tdmc_thickness <- pairwise_sma(
  model_formula = "tdmc_mean ~ thickness_mm_mean", 
  trait1 = tdmc_mean,
  trait2 = thickness_mm_mean,
  data = "log"
)

pair_sta_height <- pairwise_sma(
  model_formula = "sta_mean ~ maximum_height", 
  trait1 = sta_mean,
  trait2 = maximum_height,
  data = "log"
)

pair_sta_mh <- pairwise_sma(
  model_formula = "sta_mean ~ mass_to_height", 
  trait1 = sta_mean,
  trait2 = mass_to_height,
  data = "log"
)

pair_sta_sav <- pairwise_sma(
  model_formula = "sta_mean ~ sav_mean", 
  trait1 = sta_mean,
  trait2 = sav_mean,
  data = "log"
)

pair_sta_thickness <- pairwise_sma(
  model_formula = "sta_mean ~ thickness_mm_mean", 
  trait1 = sta_mean,
  trait2 = thickness_mm_mean,
  data = "log"
)

pair_sta_tdmc <- pairwise_sma(
  model_formula = "sta_mean ~ tdmc_mean", 
  trait1 = sta_mean,
  trait2 = tdmc_mean,
  data = "log"
)

pair_sap_height <- pairwise_sma(
  model_formula = "sap_mean ~ maximum_height", 
  trait1 = sap_mean,
  trait2 = maximum_height,
  data = "log"
)

pair_sap_sav <- pairwise_sma(
  model_formula = "sap_mean ~ sav_mean", 
  trait1 = sap_mean,
  trait2 = sav_mean,
  data = "log"
)

pair_sap_thickness <- pairwise_sma(
  model_formula = "sap_mean ~ thickness_mm_mean", 
  trait1 = sap_mean,
  trait2 = thickness_mm_mean,
  data = "log"
)

pair_sap_tdmc <- pairwise_sma(
  model_formula = "sap_mean ~ tdmc_mean", 
  trait1 = sap_mean,
  trait2 = tdmc_mean,
  data = "log"
)

pair_sap_sta <- pairwise_sma(
  model_formula = "sap_mean ~ sta_mean", 
  trait1 = sap_mean,
  trait2 = sta_mean,
  data = "log"
)     

pair_fvfm_height <- pairwise_sma(
  model_formula = "fvfm_mean ~ maximum_height", 
  trait1 = fvfm_mean,
  trait2 = maximum_height,
  data = "log"
)     

pair_fvfm_sav <- pairwise_sma(
  model_formula = "fvfm_mean ~ sav_mean", 
  trait1 = fvfm_mean,
  trait2 = sav_mean,
  data = "log"
)    

pair_fvfm_tdmc <- pairwise_sma(
  model_formula = "fvfm_mean ~ tdmc_mean", 
  trait1 = fvfm_mean,
  trait2 = tdmc_mean,
  data = "log"
)

pair_fvfm_sta <- pairwise_sma(
  model_formula = "fvfm_mean ~ sta_mean", 
  trait1 = fvfm_mean,
  trait2 = sta_mean,
  data = "log"
)

pair_fvfm_sap <- pairwise_sma(
  model_formula = "fvfm_mean ~ sap_mean", 
  trait1 = fvfm_mean,
  trait2 = sap_mean,
  data = "log"
)

pair_ar_height <- pairwise_sma(
  model_formula = "aspect_ratio_mean ~ maximum_height", 
  trait1 = aspect_ratio_mean,
  trait2 = maximum_height,
  data = "log"
)

# this also doesn't work!
# pair_ar_mh <- pairwise_sma(
#   model_formula = "aspect_ratio_mean ~ mass_to_height",
#   trait1 = aspect_ratio_mean,
#   trait2 = mass_to_height,
#   data = "log"
# )

pair_ar_sav <- pairwise_sma(
  model_formula = "aspect_ratio_mean ~ sav_mean", 
  trait1 = aspect_ratio_mean,
  trait2 = sav_mean,
  data = "log"
)

pair_ar_thickness <- pairwise_sma(
  model_formula = "aspect_ratio_mean ~ thickness_mm_mean", 
  trait1 = aspect_ratio_mean,
  trait2 = thickness_mm_mean,
  data = "log"
)

pair_ar_tdmc <- pairwise_sma(
  model_formula = "aspect_ratio_mean ~ tdmc_mean", 
  trait1 = aspect_ratio_mean,
  trait2 = tdmc_mean,
  data = "log"
)

pair_ar_fvfm <- pairwise_sma(
  model_formula = "aspect_ratio_mean ~ fvfm_mean", 
  trait1 = aspect_ratio_mean,
  trait2 = fvfm_mean,
  data = "log"
)

mh_row_scaled <- plot_grid(
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  nrow = 1
)

mh_row <- plot_grid(
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  nrow = 1
)

sav_row_scaled <- plot_grid(
  NULL,
  pair_sav_mh[[3]],
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  nrow = 1
)

sav_row <- plot_grid(
  pair_sav_height[[3]],
  pair_sav_mh[[3]],
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  nrow = 1
)

thickness_row_scaled <- plot_grid(
  NULL,
  pair_thickness_mh[[3]],
  pair_thickness_sav[[3]],
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  nrow = 1
)

thickness_row <- plot_grid(
  NULL,
  NULL,
  pair_thickness_sav[[3]],
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  nrow = 1
)

tdmc_row_scaled <- plot_grid(
  NULL,
  pair_tdmc_mh[[3]],
  NULL,
  pair_tdmc_thickness[[3]],
  NULL,
  NULL,
  NULL,
  NULL,
  nrow = 1
)

tdmc_row <- plot_grid(
  pair_tdmc_height[[3]],
  pair_tdmc_mh[[3]],
  NULL,
  pair_tdmc_thickness[[3]],
  NULL,
  NULL,
  NULL,
  NULL,
  nrow = 1
)

sta_row_scaled <- plot_grid(
  NULL,
  pair_sta_mh[[3]],
  pair_sta_sav[[3]],
  pair_sta_thickness[[3]],
  pair_sta_tdmc[[3]],
  NULL,
  NULL,
  NULL,
  nrow = 1
)

sta_row <- plot_grid(
  pair_sta_height[[3]],
  pair_sta_mh[[3]],
  pair_sta_sav[[3]],
  pair_sta_thickness[[3]],
  pair_sta_tdmc[[3]],
  NULL,
  NULL,
  NULL,
  nrow = 1
)

sap_row_scaled <- plot_grid(
  NULL,
  NULL,
  pair_sap_sav[[3]],
  pair_sap_thickness[[3]],
  pair_sap_tdmc[[3]],
  pair_sap_sta[[3]],
  NULL,
  NULL,
  nrow = 1
)

sap_row <- plot_grid(
  pair_sap_height[[3]],
  NULL,
  pair_sap_sav[[3]],
  NULL,
  pair_sap_tdmc[[3]],
  pair_sap_sta[[3]],
  NULL,
  NULL,
  nrow = 1
)

fvfm_row_scaled <- plot_grid(
  pair_fvfm_height[[3]],
  NULL,
  pair_fvfm_sav[[3]],
  NULL,
  pair_fvfm_tdmc[[3]],
  pair_fvfm_sta[[3]],
  pair_fvfm_sap[[3]],
  NULL,
  nrow = 1
)

fvfm_row <- plot_grid(
  pair_fvfm_height[[3]],
  NULL,
  pair_fvfm_sav[[3]],
  NULL,
  pair_fvfm_tdmc[[3]],
  pair_fvfm_sta[[3]],
  pair_fvfm_sap[[3]],
  NULL,
  nrow = 1
)

ar_row_scaled <- plot_grid(
  NULL,
  NULL,
  pair_ar_sav[[3]], 
  pair_ar_thickness[[3]],
  NULL,
  NULL,
  NULL,
  NULL,
  nrow = 1
)

ar_row <- plot_grid(
  pair_ar_height[[3]],
  NULL,
  pair_ar_sav[[3]], 
  pair_ar_thickness[[3]],
  pair_ar_tdmc[[3]],
  NULL,
  NULL,
  pair_ar_fvfm[[3]],
  nrow = 1
)

sma_grid_scaled <- plot_grid(mh_row_scaled,
                             sav_row_scaled,
                             thickness_row_scaled,
                             tdmc_row_scaled,
                             sta_row_scaled,
                             sap_row_scaled,
                             fvfm_row_scaled,
                             ar_row_scaled,
                             nrow = 8)

sma_grid <- plot_grid(mh_row,
                      sav_row,
                      thickness_row,
                      tdmc_row,
                      sta_row,
                      sap_row,
                      fvfm_row,
                      ar_row,
                      nrow = 8)

cowplot::save_plot(here::here(
  "figures",
  "tradeoffs",
  paste0("sma-grid_", today(), ".jpg")),
  base_height = 10,
  sma_grid_scaled)

ggsave(filename = here::here(
  "figures",
  "tradeoffs",
  paste0("sma-grid_", today(), ".jpg")),
  plot = sma_grid,
  width = 18,
  height = 18,
  units = "cm",
  dpi = 300
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -------------------- 4. Principal Components Analysis -------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section includes a Principal Components Analysis (PCA) to visualize the
# multivariate relationships between traits.

# I got a lot of help from David Zeleny's blog post: 
# https://www.davidzeleny.net/anadat-r/doku.php/en:pca_examples).


# ⟞ a. PCA ----------------------------------------------------------------

# trait by species PCA
tbs_pca <- rda(pca_mat_log)

# create a screeplot to visualize axis contributions
screeplot(tbs_pca, bstick = TRUE)

# look at the summary
summary(tbs_pca)

# proportion variance explained for downstream figure making
prop_PC1 <- "43.1%"
prop_PC2 <- "26.2%"
prop_PC3 <- "15.4%"

# ⟞ b. loadings -----------------------------------------------------------

# get loadings into data frame
loadings_df <- scores(tbs_pca, 
                      display = 'species', 
                      scaling = 0, 
                      choices = c(1, 2, 3)) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("trait") %>% 
  # arrange the data frame in order of PC1 loadings
  arrange(PC1) %>% 
  # set the factor levels so that all the traits appear in the same order
  mutate(trait = fct_inorder(trait))

# aesthetics for the loadings plots
loadings_plot_aes <- list(
  geom_col(color = "#000000", 
           fill = "orange", 
           linewidth = 0.5),
  geom_vline(xintercept = 0),
  scale_x_continuous(limits = c(-0.89, 1), 
                     breaks = seq(from = -0.75, to = 1, by = 0.25))
)

# theme elements for the loadings plots
loadings_plot_theme <- function() {
  theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 20),
          axis.title = element_blank(),
          # plot.title.position = "plot",
          plot.title = element_text(hjust = 0.5, size = 20))
}

pc1_plot <- ggplot(data = loadings_df, 
                   aes(x = PC1,
                       y = trait)) +
  loadings_plot_aes +
  loadings_plot_theme() +
  labs(title = "PC1")

pc2_plot <- ggplot(data = loadings_df, 
                   aes(x = PC2,
                       y = trait)) +
  loadings_plot_aes +
  loadings_plot_theme() +
  labs(title = "PC2")

pc3_plot <- ggplot(data = loadings_df, 
                   aes(x = PC3,
                       y = trait)) +
  loadings_plot_aes +
  loadings_plot_theme() +
  labs(title = "PC3")

loadings_plot <- pc1_plot / pc2_plot / pc3_plot

# ggsave(here::here(
#   "figures",
#   "ordination",
#   paste0("loadings_", today(), ".jpg")),
#   loadings_plot,
#   width = 14,
#   height = 16,
#   units = "cm",
#   dpi = 150
# )

# ⟞ c. PCA plots ----------------------------------------------------------

# simple biplot (compare with ggplot output to make sure it's right)
biplot(tbs_pca)

# species points
PCAscores <- scores(tbs_pca, 
                    display = "sites", 
                    choices = c(1, 2, 3)) %>% 
  as.data.frame() %>% 
  rownames_to_column("specimen_ID") %>% 
  left_join(., metadata_ind, by = "specimen_ID") %>% 
  left_join(., coarse_traits, by = "sp_code") %>% 
  mutate(scientific_name = factor(scientific_name, 
                                  levels = algae_proposal_sciname_factors))

# trait vectors
PCAvect <- scores(tbs_pca, 
                  display = "species", 
                  choices = c(1, 2, 3)) %>% 
  as.data.frame()

PCA_aesthetics <- list(
  coord_cartesian(),
    geom_vline(xintercept = c(0), color = "grey70", linetype = 2),
    geom_hline(yintercept = c(0), color = "grey70", linetype = 2),
  scale_color_manual(values = algae_colors),
  guides(colour = guide_legend(ncol = 4),
         shape = guide_legend(ncol = 4)) 
)

PCA_theme <- function() {
  theme_bw() +
    theme(legend.position = "bottom", 
          legend.box = "vertical", 
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          title = element_text(size = 20),
          legend.text = element_text(size = 10)) 
}

# plot PCA
plot_PCA_12 <- ggplot() +
  PCA_aesthetics +
  geom_point(data = PCAscores, 
             aes(x = PC1, 
                 y = PC2, 
                 color = scientific_name), 
             size = 1) +
  geom_segment(data = PCAvect, 
               aes(x = 0, 
                   y = 0, 
                   xend = PC1, 
                   yend = PC2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               linewidth = 0.5) +
  geom_text_repel(data = PCAvect, 
                  aes(x = PC1, 
                      y = PC2, 
                      label = rownames(PCAvect)), 
                  size = 8, 
                  alpha = 0.8,
                  color = "red") +
  # scale_x_continuous(limits = c(-1.75, 3.75)) +
  # scale_y_continuous(limits = c(-2, 3.25)) +
  PCA_theme() +
  labs(x = paste0("PC1 (", prop_PC1, ")"),
       y = paste0("PC2 (", prop_PC2, ")"),
       title = "PC1 and PC2",
       color = "Scientific name") 
plot_PCA_12

plot_PCA_13 <- ggplot() +
  PCA_aesthetics +
  geom_point(data = PCAscores, 
             aes(x = PC1, 
                 y = PC3, 
                 color = scientific_name), 
             size = 1) +
  geom_segment(data = PCAvect, 
               aes(x = 0, 
                   y = 0, 
                   xend = PC1, 
                   yend = PC3), 
               arrow = arrow(length = unit(0.2, "cm")), 
               linewidth = 0.5) +
  geom_text_repel(data = PCAvect, 
                  aes(x = PC1, 
                      y = PC3, 
                      label = rownames(PCAvect)), 
                  size = 8, 
                  alpha = 0.8,
                  color = "red") +
  # scale_x_continuous(limits = c(-1.75, 3.75)) +
  # scale_y_continuous(limits = c(-2, 3.25)) +
  PCA_theme() +
  labs(x = paste0("PC1 (", prop_PC1, ")"),
       y = paste0("PC3 (", prop_PC3, ")"),
       title = "PC1 and PC3",
       color = "Scientific name") 
plot_PCA_13

pca_together <- (plot_PCA_12 | plot_PCA_13) +
  plot_layout(guides = "collect", ncol = 2, widths = c(1, 1)) &
  theme(legend.position = "bottom")

pca_together

ggsave(here::here("figures", 
                  "ordination",
                  paste("PCA-log_", today(), ".jpg", sep = "")),
       pca_together,
       width = 18, height = 12, units = "cm", dpi = 300)

