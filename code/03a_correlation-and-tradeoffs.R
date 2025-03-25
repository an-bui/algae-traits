
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 0. source -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# only need to run this once per session
source(here::here("code", "01a_trait-cleaning.R"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------ 1. wrangling -----------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pca_mat <- ind_traits_filtered %>% 
  # H, T, SA, H:WW, DW:WW, H:V, SA:V, SA:DW, and SA:P
  select(specimen_ID, 
         maximum_height, 
         thickness_mm_mean, 
         # frond_area_scaled,
         height_ww,
         total_dmc, 
         height_vol
         # sav_scaled, 
         # sta_scaled,
         # sap_mean
         ) %>% 
  column_to_rownames("specimen_ID") %>% 
  drop_na() %>% 
  rename(`Height` = maximum_height,
         `Thickness` = thickness_mm_mean,
         # `Surface area` = frond_area_scaled,
         `Height:wet weight` = height_ww,
         `Dry:wet weight` = total_dmc,
         `Height:volume` = height_vol
         # `Surface area:volume` = sav_scaled,
         # `Surface area:dry weight` = sta_scaled,
         # `Surface area:perimeter` = sap_mean
         )  

pca_mat_scale <- scale(pca_mat)

pca_mat_log <- pca_mat %>% 
  # only log transforming traits that were not normally distributed
  mutate(across(where(is.numeric), 
                log))

# pca_mat_with_cn <- ind_traits %>% 
#   filter(sp_code %in% algae_proposal) %>% 
#   # weird Nienburgia?
#   filter(!(specimen_ID %in% c("20210721-BULL-023", "20210719-IVEE-009"))) %>% 
#   mutate(growth_form_num = case_when(
#     growth_form == "leathery_macrophyte" ~ 1,
#     growth_form == "corticated_macrophytes" ~ 2,
#     growth_form == "corticated_foliose" ~ 3
#   ), 
#   pigment_type_num = case_when(
#     pigment_type == "red" ~ 1,
#     pigment_type == "brown" ~ 2
#   ), 
#   longevity_num = case_when(
#     longevity == "annual" ~ 1,
#     longevity == "perennial" ~ 3,
#     longevity == "annual or perennial" ~ 2
#   )) %>% 
#   select(specimen_ID, 
#          maximum_height, total_volume, sav_mean, thickness_mm_mean, tdmc_mean, 
#          sta_mean, sav_mean, sap_mean, fvfm_mean, thickness_mm_mean, 
#          weight_percent_n_mean, weight_percent_c_mean) %>% 
#   column_to_rownames("specimen_ID") %>% 
#   drop_na()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# --------------------------- 2. plot aesthetics --------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ⟞ a. SMA plots ----------------------------------------------------------


# ⟞ b. loadings plots -----------------------------------------------------

# aesthetics for the loadings plots
loadings_plot_aes <- list(
  geom_col(color = "#000000", 
           fill = "orange", 
           linewidth = 0.5),
  geom_vline(xintercept = 0),
  scale_x_continuous(limits = c(-0.6, 0.6), 
                     breaks = seq(from = -0.6, to = 0.6, by = 0.3))
)

# theme elements for the loadings plots
loadings_plot_theme <- function() {
  theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 18),
          axis.title = element_blank(),
          # plot.title.position = "plot",
          plot.title = element_text(hjust = 0.5, size = 20))
}


# ⟞ c. PCA plots ----------------------------------------------------------

PCA_aesthetics <- list(
  # coord_cartesian(),
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2),
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2)
)

vector_colors <- list(
  scale_color_manual(values = trait_color_palette),
  scale_fill_manual(values = trait_color_palette)
)

species_colors <- list(
  scale_color_manual(values = algae_splabel_colors), 
    guides(color = guide_legend(ncol = 2))
)

PCA_theme <- function() {
  theme_bw() +
    theme(legend.position = "none", 
          plot.title.position = "plot",
          # legend.box = "vertical", 
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          title = element_text(size = 20),
          legend.text = element_text(size = 14),
          panel.grid = element_blank()) 
}

# ⟞ d. trait contributions plots ------------------------------------------

contrib_theme <- function() {
  theme_bw() +
    theme(legend.position = "none",
          plot.title.position = "plot",
          panel.grid = element_blank(),
          axis.text = element_text(size = 18),
          # axis.title = element_text(size = 20),
          plot.title = element_text(size = 22),
          axis.title = element_blank())
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ----------------------------- 3. correlations ---------------------------
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

# jpeg(here::here("figures",
#                 "correlation",
#                 paste0("corrplot_log_full-model_", today(), ".jpeg")),
#      width = 24, height = 24, units = "cm", res = 300)
# corrplot(corr = M,
#          p.mat = p$p,
#          diag = FALSE,
#          method = "circle",
#          # addgrid.col = NA,
#          type = "lower",
#          insig = "blank",
#          addCoef.col = "black",
#          col= colorRampPalette(c("#F21A00", "#FFFFFF", "#3B9AB2"))(200),
#          sig.level = 0.05)
# dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------ 4. standardized major axis regression ----------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. identify significant correlations ----------------------------------

# p-values from correlation in matrix form
p_df <- p$p

# thank you rguha for this solution: 
# https://stackoverflow.com/questions/24554642/what-is-the-best-way-to-remove-items-above-a-diagonal-in-a-data-frame-in-r

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
  pivot_longer(cols = 2:length(colnames(.)),
               names_to = "trait2",
               values_to = "corr") %>% 
  # drop any NAs in the p-values (the upper triangle of the matrix)
  drop_na(corr) %>% 
  # only include p-values that are more than 0 (the matrix has the same-same pair as 0)
  filter(corr > 0) 

log_ind_traits <- ind_traits_filtered %>% 
  mutate(across(c(maximum_height, 
                  thickness_mm_mean, 
                  # frond_area_scaled,
                  height_ww,
                  total_dmc, 
                  height_vol
                  # sav_scaled, 
                  # sta_scaled,
                  # sap_mean
                  ), 
                log)) %>% 
  mutate(sp_label = fct_relevel(sp_label, algae_splabel_factors))

pairwise_sma <- function(model_formula, trait1, trait2, data) {
  df <- if(data == "normal") {
    ind_traits_filtered
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
                     y = {{ trait2 }},
                     color = sp_label)) +
    geom_point(alpha = 0.3,
               shape = 21,
               size = 1) +
    stat_ma_line(method = "SMA",
                 color = "black",
                 fill = "lightgrey",
                 linewidth = 1) +
    geom_smooth(aes(group = sp_label),
                method = "lm",
                se = FALSE,
                linewidth = 1) +
    scale_color_manual(values = algae_splabel_colors,
                       labels = function(x) str_wrap(x, width = 40)) +
    labs(color = "Scientific name") + 
    guides(color = guide_legend(nrow = 3),
           label.position = "bottom") +
    # labs(title = model_formula) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          text = element_text(size = 26))
  
  # return all
  return(list(lmodel2_obj, sma_obj, plot))
  
  
}

pair_h_hv <- pairwise_sma(
  model_formula = "maximum_height ~ height_vol", 
  trait1 = maximum_height,
  trait2 = height_vol,
  data = "log"
)

pair_h_hv_plot <- pair_h_hv[[3]] +
  labs(x = "Height:volume ratio",
       y = "Maximum height",
       title = "(a)") +
  theme(plot.title.position = "plot")

pair_thick_height <- pairwise_sma(
  model_formula = "thickness_mm_mean ~ maximum_height", 
  trait1 = thickness_mm_mean,
  trait2 = maximum_height,
  data = "log"
)

thick_height_plot <- pair_thick_height[[3]] +
  labs(x = "Height",
       y = "Thickness",
       title = "(c)") +
  theme(plot.title.position = "plot")

pair_t_h_ww <- pairwise_sma(
  model_formula = "thickness_mm_mean ~ height_ww", 
  trait1 = thickness_mm_mean,
  trait2 = height_ww,
  data = "log"
)

pair_t_h_ww <- pair_t_h_ww[[3]] +
  labs(x = "Height:wet weight",
       y = "Thickness",
       title = "(b)") +
  theme(plot.title.position = "plot")

pair_dmc_height <- pairwise_sma(
  model_formula = "total_dmc ~ maximum_height", 
  trait1 = total_dmc,
  trait2 = maximum_height,
  data = "log"
)

dmc_height_plot <- pair_dmc_height[[3]] +
  labs(x = "Height",
       y = "Dry:wet weight",
       title = "(d)") +
  theme(plot.title.position = "plot")

plot_legend <- pair_dmc_height[[3]] +
  labs(color = "Scientific name") + 
  guides(color = guide_legend(nrow = 3),
         label.position = "bottom")

plot_legend_test <- cowplot::get_plot_component(plot_legend, "guide-box-right", return_all = TRUE)


sma_together <- (pair_h_hv_plot | thick_height_plot) / (pair_t_h_ww | dmc_height_plot) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  
# ggsave(here::here("figures",
#                   "tradeoffs",
#                   paste0("sma_supplement_", today(), ".jpg")),
#        width = 16,
#        height = 18,
#        units = "cm",
#        dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -------------------- 5. Principal Components Analysis -------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section includes a Principal Components Analysis (PCA) to visualize the
# multivariate relationships between traits.

# I got a lot of help from David Zeleny's blog post: 
# https://www.davidzeleny.net/anadat-r/doku.php/en:pca_examples).


# ⟞ a. PCA ----------------------------------------------------------------

# trait by species PCA
pca_full <- rda(pca_mat_log, scale = TRUE)

# create a screeplot to visualize axis contributions
screeplot(pca_full, bstick = TRUE)

# look at the summary
summary(pca_full)
prcomp(pca_mat_log, center = TRUE, scale. = TRUE) %>% summary()

# proportion variance explained for downstream figure making
prop_PC1_full <- "53.2%"
prop_PC2_full <- "28.4%"

trait_dist <- vegdist(pca_mat_log, 
                      method = "euclidean")

sp_permanova <- adonis2(trait_dist ~ sp_code, 
                        data = ind_traits_filtered)

# sp_permanova %>%
#   tidy() %>%
#   mutate(across(SumOfSqs:statistic, ~ round(.x, digits = 2))) %>%
#   flextable() %>%
#   autofit() %>%
#   fit_to_width(8) %>%
#   font(fontname = "Times New Roman",
#        part = "all") %>%
#   # change the column names
#   set_header_labels("term" = "Term",
#                     "df" = "Degrees of freedom",
#                     "SumOfSqs" = "Sum of squares",
#                     "R2" = "R\U00B2",
#                     "statistic" = "F-statistic",
#                     "p.value" = "p-value") %>%
#   save_as_docx(path = here("tables", "PERMANOVA", paste0("full-trait-ANOVA_", today(), ".docx")))
# species are different from each other

set.seed(1)

rvam_pairwise_log <- pairwise.perm.manova(resp = trait_dist,
                                          fact = ind_traits_filtered$sp_code,
                                          p.method = "BH",
                                          `F` = TRUE)

rvam_pairwise_log
# BO and CO are not different from each other = articulated corallines
# DP and BF are not different = foliose ish things
# PTCA and LAFA are not different = stipate browns
# CYOS: bushy stipate brown
# CC: short leathery red
# R: bushy short red
# DP: bushy short brown

full_pairwise_matrix <- (rvam_pairwise_log$p.value < 0.05)

# rvam_pairwise_log$p.value %>%
#   as.data.frame() %>%
#   rownames_to_column("sp_code") %>%
#   mutate(across(c(CO:PTCA), ~ case_when(
#     between(.x, 0, 0.001) ~ "<0.001",
#     between(.x, 0.001, 0.01) ~ as.character(round(.x, digits = 3)),
#     between(.x, 0.01, 1) ~ as.character(round(.x, digits = 2))
#   ))) %>%
#   mutate(across(everything(), ~ replace_na(.x, "-"))) %>%
#   flextable() %>%
#   autofit() %>%
#   fit_to_width(8) %>%
#   font(fontname = "Times New Roman",
#        part = "all") %>%
#   set_header_labels("sp_code" = "") %>%
#   save_as_docx(path = here("tables", "PERMANOVA", paste0("full-trait-ANOVA_pairwise-comparisons_", today(), ".docx")))

anova(betadisper(d = trait_dist,
                 group = ind_traits_filtered$sp_code)) %>% 
  tidy() %>% 
  mutate(across(where(is.numeric), ~round(., digits = 2))) %>% 
  flextable() %>% 
  set_header_labels("term" = "Term",
                    "df" = "Degrees of freedom",
                    "sumsq" = "Sum of squares",
                    "meansq" = "Mean squares",
                    "statistic" = "F-statistic",
                    "p.value" = "p-value") %>% 
  autofit() %>% 
  fit_to_width(5, unit = "in") %>% 
  font(fontname = "Times New Roman",
       part = "all") %>% 
  save_as_docx(path = here("tables",
                    "PERMANOVA",
                    paste0("full-trait_dispersions_", today(), ".docx")))
  
# no difference in dispersions

# ⟞ b. loadings -----------------------------------------------------------

# get loadings into data frame
loadings_df_full <- scores(pca_full, 
                           display = 'species', 
                           scaling = 0, 
                           choices = c(1, 2)) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("trait") %>% 
  # arrange the data frame in order of PC1 loadings
  arrange(PC1) %>% 
  # set the factor levels so that all the traits appear in the same order
  mutate(trait = fct_inorder(trait))

pc1_plot_full <- ggplot(data = loadings_df_full, 
                   aes(x = PC1,
                       y = trait)) +
  loadings_plot_aes +
  loadings_plot_theme() +
  labs(title = "PC1")

pc2_plot_full <- ggplot(data = loadings_df_full, 
                   aes(x = PC2,
                       y = trait)) +
  loadings_plot_aes +
  loadings_plot_theme() +
  labs(title = "PC2")

loadings_plot_full <- pc1_plot_full / pc2_plot_full 

# ggsave(here::here(
#   "figures",
#   "ordination",
#   paste0("loadings_scale_full-model_", today(), ".jpg")),
#   loadings_plot_full,
#   width = 12,
#   height = 14,
#   units = "cm",
#   dpi = 300
# )

# ⟞ c. biplots ------------------------------------------------------------

# simple biplot (compare with ggplot output to make sure it's right)
biplot(pca_full)

# species points
PCAscores_full <- scores(pca_full, 
                    display = "sites", 
                    choices = c(1, 2)) %>% 
  as.data.frame() %>% 
  rownames_to_column("specimen_ID") %>% 
  left_join(., metadata_ind, by = "specimen_ID") %>% 
  left_join(., coarse_traits, by = "sp_code") %>% 
  mutate(scientific_name = factor(scientific_name, 
                                  levels = algae_factors)) %>% 
  left_join(., fvfm_ind, by = "specimen_ID") %>% 
  left_join(., enframe(algae_spcode_full_names), by = c("sp_code" = "name")) %>% 
  rename("splabel" = "value") %>% 
  relocate(splabel, .after = scientific_name) %>% 
  mutate(splabel = factor(splabel,
                          levels = algae_splabel_factors))


# ⟞ ⟞ i. trait vectors ----------------------------------------------------

# trait vectors
PCAvect_full <- scores(pca_full, 
                       display = "species", 
                       choices = c(1, 2)) %>% 
  as.data.frame()

# plot PCA
plot_PCA_12_vectors <- ggplot() +
  PCA_aesthetics +
  vector_colors + 
  geom_point(data = PCAscores_full, 
             aes(x = PC1, 
                 y = PC2#, 
                 # color = scientific_name, 
                 # shape = scientific_name,
                 # size = fvfm_mean
                 ) ,
             alpha = 0.3,
             size = 1,
             shape = 21,
             color = "darkgrey"
             ) +
  geom_segment(data = PCAvect_full, 
               aes(x = 0, 
                   y = 0, 
                   xend = PC1, 
                   yend = PC2,
                   color = rownames(PCAvect_full)), 
               arrow = arrow(length = unit(0.2, "cm")), 
               linewidth = 0.5) +
  geom_label_repel(data = PCAvect_full, 
                  aes(x = PC1, 
                      y = PC2, 
                      label = rownames(PCAvect_full),
                      fill = rownames(PCAvect_full)), 
                  size = 6, 
                  alpha = 0.8,
                  seed = 666,
                  color = "black") +
  scale_x_continuous(limits = c(-3, 3)) +
  scale_y_continuous(limits = c(-3, 3)) +
  PCA_theme() +
  labs(x = paste0("PC1 (", prop_PC1_full, ")"),
       y = paste0("PC2 (", prop_PC2_full, ")"),
       title = "(a) Trait vectors",
       # subtitle = "BO, CC, CO, BF, DP, LAFA, PTCA, R, CYOS", 
       color = "Scientific name") 
plot_PCA_12_vectors

# ggsave(here::here("figures",
#                   "ordination",
#                   paste("PCA-log_scale_full-model_vectors_", today(), ".jpg", sep = "")),
#        plot_PCA_12_vectors,
#        width = 8, height = 8, units = "cm", dpi = 300)


# ⟞ ⟞ i. species points ---------------------------------------------------

plot_PCA_12_species <- PCAscores_full %>% 
  ggplot(aes(x = PC1, 
             y = PC2, 
             color = splabel, 
             # shape = scientific_name,
             # size = fvfm_mean
  ) ) +
  PCA_aesthetics +
  species_colors + 
  geom_point(
    shape = 21,
    alpha = 0.3,
    size = 1
  ) +
  stat_ellipse(aes(color = splabel),
               level = 0.5) +
  scale_x_continuous(limits = c(-1.05, 1.05)) +
  scale_y_continuous(limits = c(-1.05, 1.05)) +
  PCA_theme() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.51, 0.89),
        legend.background = element_blank(),
        legend.spacing.y = unit(0.01, "cm"),
        legend.key.spacing.y = unit(0.01, "cm"),
        legend.key.height = unit(0.25, "cm"),
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_blank()) +
  labs(color = "Scientific name",
       title = "(b) Species",
       x = paste0("PC1 (", prop_PC1_full, ")"),
       y = paste0("PC2 (", prop_PC2_full, ")"),
  )

plot_PCA_12_species

PCA_vectors_species <- plot_PCA_12_vectors | plot_PCA_12_species

# ggsave(here::here("figures",
#                   "ordination",
#                   paste("PCA_full-model_vectors-and-species_", today(), ".jpg", sep = "")),
#        PCA_vectors_species,
#        width = 16, height = 8, units = "cm", dpi = 300)

# ⟞ d. axis contributions -------------------------------------------------

# got calculation from http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
varcoord_full <- PCAvect_full %>% 
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
  mutate(contrib = (values/component_total)*100) %>% 
  left_join(., enframe(trait_abbreviations), by = c("trait" = "value")) %>% 
  rename("abbrev" = "name")

# The red dashed line on the graph above indicates the expected average contribution. If the contribution of the variables were uniform, the expected value would be 1/length(variables) = 1/9 = 11.1%
expected_average_full <- (1/length(unique(varcoord_full$trait)))*100

contrib_aesthetics_full <- list(
  geom_col(color = "black"),
  geom_hline(yintercept = expected_average_full,
             color = "#96a39e",
             linetype = 2),
  scale_fill_manual(values = trait_color_palette),
  scale_y_continuous(expand = expansion(c(0, 0.05)),
                     position = "left")
)

pc1_contrib_full <- varcoord_full %>% 
  filter(axis == "PC1") %>% 
  ggplot(aes(x = reorder(trait, -contrib),
             y = contrib, 
             fill = trait)) +
  contrib_aesthetics_full +
  contrib_theme() +
  labs(title = "(c) Trait % contributions to PC1")

pc1_contrib_full

pc2_contrib_full <- varcoord_full %>% 
  filter(axis == "PC2") %>% 
  ggplot(aes(x = reorder(trait, -contrib),
             y = contrib, 
             fill = trait)) +
  contrib_aesthetics_full +
  contrib_theme() +
  labs(title = "(d) Trait % contributions to PC2")

pc2_contrib_full

# scaled layout
contrib_together_full <- 
  (
    plot_PCA_12_vectors | plot_PCA_12_species
  ) / (
    pc1_contrib_full / pc2_contrib_full
  ) +
  plot_layout(axis_titles = "collect")

# contrib_together_full <- (free(plot_PCA_12_full) / 
#                             (free(pc1_contrib_full) | free(pc2_contrib_full))) + 
#   plot_layout(heights = c(1, 0.8))

# contrib_together_full <- ((free(plot_PCA_12_full) / (free(pc1_contrib_full)) | (free(pc2_contrib_full) / plot_spacer()))) + 
#   plot_layout(heights = c(1, 0.5, 0.5))

# ggsave(here::here(
#   "figures",
#   "ordination",
#   paste0("contributions_scale_full-model_", today(), ".jpg")),
#   contrib_together_full,
#   width = 16,
#   height = 16,
#   units = "cm",
#   dpi = 300
# )





