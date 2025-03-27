
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 0. source -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# only need to run this once per session
source(here::here("code", "01a_trait-cleaning.R"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------ 1. wrangling -----------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# old matrix for testing
# pca_mat <- ind_traits_filtered %>% 
#   # H, T, SA, H:WW, DW:WW, H:V, SA:V, SA:DW, and SA:P
#   select(specimen_ID, 
#          maximum_height, 
#          thickness_mm_mean, 
#          # frond_area_scaled,
#          height_ww,
#          total_dmc, 
#          height_vol
#          # sav_scaled, 
#          # sta_scaled,
#          # sap_mean
#   ) %>% 
#   column_to_rownames("specimen_ID") %>% 
#   drop_na() %>% 
#   rename(`Height` = maximum_height,
#          `Thickness` = thickness_mm_mean,
#          #`Surface area` = frond_area_scaled,
#          `Height:wet weight` = height_ww,
#          `Dry:wet weight` = total_dmc,
#          `Height:volume` = height_vol
#          # `Surface area:volume` = sav_scaled,
#          # `Surface area:dry weight` = sta_scaled,
#          # `Surface area:perimeter` = sap_mean
#   )  
# 
# pca_mat_log <- pca_mat %>% 
#   # only log transforming traits that were not normally distributed
#   mutate(across(where(is.numeric), 
#                 log))

pca_mat_full <- ind_traits_filtered_full %>% 
  # H, T, SA, H:WW, DW:WW, H:V, SA:V, SA:DW, and SA:P
  select(specimen_ID, 
         maximum_height, 
         thickness_mm_mean, 
         frond_area_scaled,
         height_ww,
         total_dmc, 
         height_vol,
         sav_scaled, 
         sta_scaled,
         sap_mean
         ) %>% 
  column_to_rownames("specimen_ID") %>% 
  drop_na() %>% 
  rename(`Height` = maximum_height,
         `Thickness` = thickness_mm_mean,
         `Surface area` = frond_area_scaled,
         `Height:wet weight` = height_ww,
         `Dry:wet weight` = total_dmc,
         `Height:volume` = height_vol,
         `Surface area:volume` = sav_scaled,
         `Surface area:dry weight` = sta_scaled,
         `Surface area:perimeter` = sap_mean
         )  

pca_mat_full_log <- pca_mat_full %>% 
  # only log transforming traits that were not normally distributed
  mutate(across(where(is.numeric), 
                log))

pca_mat_sub <- pca_mat_full %>% 
  select(`Height`,
         `Thickness`,
         `Height:wet weight`,
         `Dry:wet weight`,
         `Height:volume`
  )

pca_mat_sub_log <- pca_mat_full_log %>% 
  select(`Height`,
         `Thickness`,
         `Height:wet weight`,
         `Dry:wet weight`,
         `Height:volume`
  )

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

## a. loadings plots ------------------------------------------------------

# aesthetics for the loadings plots
loadings_plot_aes <- list(
  geom_col(color = "#000000", 
           fill = "orange", 
           linewidth = 0.5),
  geom_vline(xintercept = 0),
  scale_x_continuous(limits = c(-1, 1), 
                     breaks = seq(from = -1, to = 1, by = 0.5))
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


## b. PCA plots -----------------------------------------------------------

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

## c. trait contributions plots -------------------------------------------

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

## a. subset traits -------------------------------------------------------

# create a correlation matrix based on pearson correlation
corr_mat_sub <- pca_mat_sub_log %>% 
  cor(method = "pearson")

# extract p-values
p_sub <- cor.mtest(pca_mat_sub_log,
                   conf.level = 0.95,
                   method = "pearson")

# create correlation plot
corr_traits_sub <- corrplot(corr = corr_mat_sub, 
                            p.mat = p_sub$p,
                            diag = FALSE,
                            method = "circle", 
                            # addgrid.col = NA,
                            type = "lower", 
                            insig = "blank", 
                            addCoef.col = "black", 
                            col= colorRampPalette(c("#F21A00", 
                                                    "#FFFFFF", 
                                                    "#3B9AB2"))(200), 
                            sig.level = 0.05)

# save figure
# jpeg(here::here("figures",
#                 "correlation",
#                 paste0("corrplot_log_sub-trait-set_", today(), ".jpeg")),
#      width = 24, height = 24, units = "cm", res = 300)
# corrplot(corr = corr_mat_sub, 
#          p.mat = p_sub$p,
#          diag = FALSE,
#          method = "circle", 
#          # addgrid.col = NA,
#          type = "lower", 
#          insig = "blank", 
#          addCoef.col = "black", 
#          col= colorRampPalette(c("#F21A00", 
#                                  "#FFFFFF", 
#                                  "#3B9AB2"))(200), 
#          sig.level = 0.05)
# dev.off()

## b. all traits ----------------------------------------------------------

# create a correlation matrix based on pearson correlation
corr_mat_full <- pca_mat_full_log %>% 
  cor(method = "pearson")

# extract p-values
p_full <- cor.mtest(pca_mat_full_log,
                    conf.level = 0.95,
                    method = "pearson")

# create correlation plot
corr_traits_full <- corrplot(corr = corr_mat_full, 
                             p.mat = p_full$p,
                             diag = FALSE,
                             method = "circle", 
                             # addgrid.col = NA,
                             type = "lower", 
                             insig = "blank", 
                             addCoef.col = "black", 
                             col= colorRampPalette(c("#F21A00", 
                                                     "#FFFFFF", 
                                                     "#3B9AB2"))(200), 
                             sig.level = 0.05)


# save figure
# jpeg(here::here("figures",
#                 "correlation",
#                 paste0("corrplot_log_full-trait-set_", today(), ".jpeg")),
#      width = 24, height = 24, units = "cm", res = 300)
# corrplot(corr = corr_mat_full, 
#          p.mat = p_full$p,
#          diag = FALSE,
#          method = "circle", 
#          # addgrid.col = NA,
#          type = "lower", 
#          insig = "blank", 
#          addCoef.col = "black", 
#          col= colorRampPalette(c("#F21A00", 
#                                  "#FFFFFF", 
#                                  "#3B9AB2"))(200), 
#          sig.level = 0.05)
# dev.off()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------ 4. standardized major axis regression ----------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# This section includes code to conduct standardized major axis regression. The
# output at the end is a figure of the predictions from the SMA and a table
# of p-values indicating significant relationships from the SMA. I created the
# p-value diagonal display manually using a solution on Stack Overflow posted
# by rauha: https://stackoverflow.com/questions/24554642/what-is-the-best-way-to-remove-items-above-a-diagonal-in-a-data-frame-in-r
# This step is important because one of the assumptions of SMA is that there is
# a significant correlation between the variables in question. So, this section
# builds on the correlation section and uses the `p_full` and `p_sub` objects
# to determine what (if any) significant correlative relationships could then
# be tested using SMA.

## a. identify significant correlations -----------------------------------

# p-values from correlation in matrix form
p_df <- p_full$p

# get rid of the upper triangle
p_df[upper.tri(p_df)] <- NA

# this creates a data frame of traits in pairs (column names: trait1 and 
# trait2), and their correlation p-value (column name: corr). Any significant 
# correlation could be tested with SMA.
pairwise_comparisons <- p_df %>%
  # turn into data frame
  as.data.frame() %>% 
  # make the rownames into "trait1" or the first trait in the pairwise comparison
  rownames_to_column("trait1") %>%
  # when p < 0.05, keep the p, and when it's not then replace with NA
  mutate(across(where(is.numeric),
                ~ case_when(. < 0.05 ~ ., TRUE ~ NA))) %>%
  # make the data frame longer and turn the column names into "trait2" 
  # or the second trait in the comparison
  pivot_longer(cols = 2:length(colnames(.)),
               names_to = "trait2",
               values_to = "corr") %>% 
  # drop any NAs in the p-values (the upper triangle of the matrix)
  drop_na(corr) %>% 
  # only include p-values that are more than 0 
  # (the matrix has the same-same pair as 0)
  filter(corr > 0) 

log_ind_traits <- ind_traits_filtered_full %>% 
  mutate(across(c(maximum_height, 
                  thickness_mm_mean, 
                  frond_area_scaled,
                  height_ww,
                  total_dmc, 
                  height_vol,
                  sav_scaled, 
                  sta_scaled,
                  sap_mean
                  ), 
                log)) %>% 
  mutate(sp_label = fct_relevel(sp_label, algae_splabel_factors))


## b. function for generating figures and analysis ------------------------

# This is a function to do the SMA and create a figure. The output is a list: 
# 1) SMA output from `lmodel2()`
# 2) SMA output from `sma()`
# 3) figure of SMA line

pairwise_sma <- function(model_formula, trait1, trait2, data) {
  df <- if(data == "normal") {
    ind_traits_filtered_full
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

## c. SMA analysis and figures --------------------------------------------

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

pair_sav_sadw <- pairwise_sma(
  model_formula = "sav_scaled ~ sta_scaled", 
  trait1 = sav_scaled,
  trait2 = sta_scaled,
  data = "log"
)

sav_sadw_plot <- pair_sav_sadw[[3]] +
  labs(x = "Surface area:volume",
       y = "Surface area:dry weight",
       title = "(a)") +
  theme(plot.title.position = "plot")

pair_hww_sadw <- pairwise_sma(
  model_formula = "height_ww ~ sta_scaled", 
  trait1 = height_ww,
  trait2 = sta_scaled,
  data = "log"
)

hww_sadw_plot <- pair_hww_sadw[[3]] +
  labs(x = "Surface area:volume",
       y = "Height:wet weight",
       title = "(b)") +
  theme(plot.title.position = "plot")


## d. putting figures together and saving ---------------------------------

sma_together <- (pair_h_hv_plot | thick_height_plot) / (pair_t_h_ww | dmc_height_plot) + 
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')

supp_sma <- (sav_sadw_plot | hww_sadw_plot) +
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')
  
# ggsave(here::here("figures",
#                   "tradeoffs",
#                   paste0("sma_app1_", today(), ".jpg")),
#        sma_together,
#        width = 16,
#        height = 18,
#        units = "cm",
#        dpi = 300)
# 
# ggsave(here::here("figures",
#                   "tradeoffs",
#                   paste0("sma_app2_", today(), ".jpg")),
#        supp_sma,
#        width = 16,
#        height = 10,
#        units = "cm",
#        dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -------------------- 5. Principal Components Analysis -------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section includes a Principal Components Analysis (PCA) to visualize the
# multivariate relationships between traits. Section a handles the subset traits
# which appear in the main text of the paper, and section b handles the full 
# set of traits.

# I got a lot of help from David Zeleny's blog post: 
# https://www.davidzeleny.net/anadat-r/doku.php/en:pca_examples).

# Additionally, the calculations for the average (expected) contributions of 
# each trait to the PCA axes is from here: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
# I couldn't figure out where it came from, so hopefully this is helpful to 
# someone out there!

## a. subset --------------------------------------------------------------

### i. PCA and PERMANOVA --------------------------------------------------

# trait by species PCA
pca_sub <- rda(pca_mat_sub_log, scale = TRUE)

# create a screeplot to visualize axis contributions
screeplot(pca_sub, bstick = TRUE)

# look at the summary
summary(pca_sub)

# proportion variance explained for downstream figure making
prop_PC1_sub <- "53.2%"
prop_PC2_sub <- "28.4%"

trait_dist_sub <- vegdist(pca_mat_sub_log, 
                      method = "euclidean")

sp_permanova <- adonis2(trait_dist ~ sp_code, 
                        data = ind_traits_filtered_sub)

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
#   save_as_docx(path = here("tables",
#                            "PERMANOVA",
#                            paste0("subset_full-trait-ANOVA_", today(), ".docx")))
# species are different from each other

set.seed(1)

pairwise_sub <- pairwise.perm.manova(resp = trait_dist_sub,
                                          fact = ind_traits_filtered_sub$sp_code,
                                          p.method = "BH",
                                          `F` = TRUE)

pairwise_sub
# BO and CO are not different from each other = articulated corallines
# DP and BF are not different = foliose ish things
# PTCA and LAFA are not different = stipate browns
# CYOS: bushy stipate brown
# CC: short leathery red
# R: bushy short red
# DP: bushy short brown

# creating this object for downstream trait selection stuff
pairwise_sub_matrix <- (pairwise_sub$p.value < 0.05)

# pairwise_sub$p.value %>%
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
#   save_as_docx(path = here("tables",
#                            "PERMANOVA",
#                            paste0(
#                              "subset_full-trait-ANOVA_pairwise-comparisons_", 
#                              today(), 
#                              ".docx")))

# anova(betadisper(d = trait_dist_sub,
#                  group = ind_traits_filtered_sub$sp_code)) %>%
#   tidy() %>%
#   mutate(across(where(is.numeric), ~round(., digits = 2))) %>%
#   flextable() %>%
#   set_header_labels("term" = "Term",
#                     "df" = "Degrees of freedom",
#                     "sumsq" = "Sum of squares",
#                     "meansq" = "Mean squares",
#                     "statistic" = "F-statistic",
#                     "p.value" = "p-value") %>%
#   autofit() %>%
#   fit_to_width(5, unit = "in") %>%
#   font(fontname = "Times New Roman",
#        part = "all") %>%
#   save_as_docx(path = here("tables",
#                     "PERMANOVA",
#                     paste0("subset_full-trait_dispersions_", today(), ".docx")))
# no difference in dispersions

### ii. loadings ----------------------------------------------------------

# get loadings into data frame
loadings_df_sub <- scores(pca_sub, 
                           display = 'species', 
                           scaling = 0, 
                           choices = c(1, 2)) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column("trait") %>% 
  # arrange the data frame in order of PC1 loadings
  arrange(PC1) %>% 
  # set the factor levels so that all the traits appear in the same order
  mutate(trait = fct_inorder(trait))

pc1_plot_sub <- ggplot(data = loadings_df_sub, 
                        aes(x = PC1,
                            y = trait)) +
  loadings_plot_aes +
  loadings_plot_theme() +
  labs(title = "PC1")

pc2_plot_sub <- ggplot(data = loadings_df_sub, 
                        aes(x = PC2,
                            y = trait)) +
  loadings_plot_aes +
  loadings_plot_theme() +
  labs(title = "PC2")

loadings_plot_sub <- pc1_plot_sub / pc2_plot_sub

# ggsave(here::here(
#   "figures",
#   "ordination",
#   paste0("loadings_scale_subset_full-model_", today(), ".jpg")),
#   loadings_plot_sub,
#   width = 12,
#   height = 14,
#   units = "cm",
#   dpi = 300
# )

### iii. biplots ----------------------------------------------------------

# simple biplot (compare with ggplot output to make sure it's right)
biplot(pca_sub)

# species points
PCAscores_sub <- scores(pca_sub, 
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

#### 1. trait vectors -----------------------------------------------------

# trait vectors
PCAvect_sub <- scores(pca_sub, 
                       display = "species", 
                       choices = c(1, 2)) %>% 
  as.data.frame()

# plot PCA
plot_PCA_12_vectors_sub <- ggplot() +
  PCA_aesthetics +
  vector_colors + 
  geom_point(data = PCAscores_sub, 
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
  geom_segment(data = PCAvect_sub, 
               aes(x = 0, 
                   y = 0, 
                   xend = PC1, 
                   yend = PC2,
                   color = rownames(PCAvect_sub)), 
               arrow = arrow(length = unit(0.2, "cm")), 
               linewidth = 0.5) +
  geom_label_repel(data = PCAvect_sub, 
                   aes(x = PC1, 
                       y = PC2, 
                       label = rownames(PCAvect_sub),
                       fill = rownames(PCAvect_sub)), 
                   size = 6, 
                   alpha = 0.8,
                   seed = 666,
                   color = "black") +
  scale_x_continuous(limits = c(-3, 3)) +
  scale_y_continuous(limits = c(-3, 3)) +
  PCA_theme() +
  labs(x = paste0("PC1 (", prop_PC1_sub, ")"),
       y = paste0("PC2 (", prop_PC2_sub, ")"),
       title = "(a) Trait vectors",
       # subtitle = "BO, CC, CO, BF, DP, LAFA, PTCA, R, CYOS", 
       color = "Scientific name") 

#### 2. species points ----------------------------------------------------

plot_PCA_12_species_sub <- PCAscores_sub %>% 
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
       x = paste0("PC1 (", prop_PC1_sub, ")"),
       y = paste0("PC2 (", prop_PC2_sub, ")"),
  )

PCA_vectors_species_sub <- plot_PCA_12_vectors_sub | plot_PCA_12_species_sub

### iv. axis contributions ------------------------------------------------

varcoord_sub <- PCAvect_sub %>% 
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
expected_average_sub <- (1/length(unique(varcoord_sub$trait)))*100

contrib_aesthetics_sub <- list(
  geom_col(color = "black"),
  geom_hline(yintercept = expected_average_sub,
             color = "#96a39e",
             linetype = 2),
  scale_fill_manual(values = trait_color_palette),
  scale_y_continuous(expand = expansion(c(0, 0.05)),
                     position = "left")
)

pc1_contrib_sub <- varcoord_sub %>% 
  filter(axis == "PC1") %>% 
  ggplot(aes(x = reorder(trait, -contrib),
             y = contrib, 
             fill = trait)) +
  contrib_aesthetics_sub +
  contrib_theme() +
  labs(title = "(c) Trait % contributions to PC1")

pc2_contrib_sub <- varcoord_sub %>% 
  filter(axis == "PC2") %>% 
  ggplot(aes(x = reorder(trait, -contrib),
             y = contrib, 
             fill = trait)) +
  contrib_aesthetics_sub +
  contrib_theme() +
  labs(title = "(d) Trait % contributions to PC2")

### v. compiling figures and saving ---------------------------------------

# scaled layout
contrib_together_sub <- 
  (
    plot_PCA_12_vectors_sub | plot_PCA_12_species_sub
  ) / (
    pc1_contrib_sub / pc2_contrib_sub
  ) +
  plot_layout(axis_titles = "collect")

# ggsave(here::here(
#   "figures",
#   "ordination",
#   paste0("contributions_subset-model_", today(), ".jpg")),
#   contrib_together_sub,
#   width = 16,
#   height = 16,
#   units = "cm",
#   dpi = 300
# )

## a. full set ------------------------------------------------------------
### i. PCA and PERMANOVA --------------------------------------------------
### ii. loadings ----------------------------------------------------------
### iii. biplots ----------------------------------------------------------
#### 1. trait vectors -----------------------------------------------------
#### 2. species points ----------------------------------------------------
### iv. axis contributions ------------------------------------------------
### v. compiling figures and saving ---------------------------------------

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





