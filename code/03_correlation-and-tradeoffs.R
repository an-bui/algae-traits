
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 0. source -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# only need to run this once per session
source(here::here("code", "01a_trait-cleaning.R"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------ 1. wrangling -----------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pca_mat <- ind_traits_filtered %>% 
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
  # select(!(`Fv/Fm`))

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
  coord_cartesian(),
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2),
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2),
  scale_color_manual(values = trait_color_palette),
  scale_fill_manual(values = trait_color_palette),
  guides(colour = guide_legend(ncol = 4),
         shape = guide_legend(ncol = 4)) 
)

PCA_theme <- function() {
  theme_bw() +
    theme(legend.position = "none", 
          plot.title.position = "plot",
          # legend.box = "vertical", 
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          title = element_text(size = 20),
          legend.text = element_text(size = 10)) 
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


# ⟞ a. full model ---------------------------------------------------------

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
                paste0("corrplot_log_full-model_", today(), ".jpeg")),
     width = 24, height = 24, units = "cm", res = 300)
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

# ⟞ b. reduced model ------------------------------------------------------

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
# ------------------ 4. standardized major axis regression ----------------
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
  pivot_longer(cols = 2:length(colnames(.)),
               names_to = "trait2",
               values_to = "corr") %>%
  # drop any NAs in the p-values (the upper triangle of the matrix)
  drop_na(corr) %>% 
  # only include p-values that are more than 0 (the matrix has the same-same pair as 0)
  filter(corr > 0) %>% 
  mutate(trait1_column = case_match(
    trait1,
    "Height" ~ "maximum_height",
    "Dry weight:height" ~ "mass_to_height",
    "Dry weight" ~ "total_dry",
    "Dry:wet weight" ~ "total_dmc",
    "Wet weight" ~ "total_wet",
    "Volume" ~ "total_volume",
    "Surface area:volume" ~ "sav_scaled",
    "Surface area" ~ "frond_area_scaled",
    "Thickness" ~ "thickness_mm_mean",
    "Surface area:dry weight" ~ "sta_scaled",
    "Surface area:perimeter" ~ "sap_mean",
    "Perimeter" ~ "frond_peri_scaled"
  )) %>% 
  mutate(trait2_column = case_match(
    trait2,
    "Height" ~ "maximum_height",
    "Dry weight:height" ~ "mass_to_height",
    "Dry weight" ~ "total_dry",
    "Dry:wet weight" ~ "total_dmc",
    "Wet weight" ~ "total_wet",
    "Volume" ~ "total_volume",
    "Surface area:volume" ~ "sav_scaled",
    "Surface area" ~ "frond_area_scaled",
    "Thickness" ~ "thickness_mm_mean",
    "Surface area:dry weight" ~ "sta_scaled",
    "Surface area:perimeter" ~ "sap_mean"
  )) %>% 
  mutate(formula = paste(trait1_column, "~", trait2_column, sep = " ")) 

log_ind_traits <- ind_traits_filtered %>% 
  mutate(across(c(maximum_height, mass_to_height, total_dry, 
                  total_dmc, total_wet, total_volume, 
                  sav_scaled, frond_area_scaled,
                  thickness_mm_mean,
                  sta_scaled, sap_mean, frond_peri_scaled), 
                log))

scaled_ind_traits <- ind_traits %>% 
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
    geom_point(alpha = 0.3,
               shape = 21,
               color = "#7e97c0") +
    stat_ma_line(method = "SMA",
                 color = "#295396",
                 fill = "#a9bad5",
                 linewidth = 0.75) +
    # labs(title = model_formula) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          text = element_text(size = 12))
  
  # return all
  return(list(lmodel2_obj, sma_obj, plot))
  
  
}

pair_sta_sav <- pairwise_sma(
  model_formula = "sta_scaled ~ sav_scaled", 
  trait1 = sta_scaled,
  trait2 = sav_scaled,
  data = "log"
)

sta_sav_plot <- pair_sta_sav[[3]] +
  labs(x = "Surface area:dry weight",
       y = "Surface area:volume ratio",
       title = "(a)") +
  theme(plot.title.position = "plot")

pair_sta_mh <- pairwise_sma(
  model_formula = "sta_scaled ~ mass_to_height", 
  trait1 = sta_scaled,
  trait2 = mass_to_height,
  data = "log"
)

sta_mh_plot <- pair_sta_mh[[3]] +
  labs(x = "Surface area:dry weight",
       y = "Mass:height",
       title = "(b)") +
  theme(plot.title.position = "plot")

pair_dmc_height <- pairwise_sma(
  model_formula = "total_dmc ~ maximum_height", 
  trait1 = total_dmc,
  trait2 = maximum_height,
  data = "log"
)

dmc_height_plot <- pair_dmc_height[[3]] +
  labs(x = "Dry:wet weight",
       y = "Height",
       title = "(c)") +
  theme(plot.title.position = "plot")

sma_together <- sta_sav_plot + sta_mh_plot + dmc_height_plot

ggsave(here::here("figures",
                  "tradeoffs",
                  paste0("sma_main-text_", today(), ".jpg")),
       width = 10,
       height = 4,
       units = "cm",
       dpi = 300)

# pair_sav_height <- pairwise_sma(
#   model_formula = "sav_mean ~ maximum_height", 
#   trait1 = sav_mean,
#   trait2 = maximum_height,
#   data = "log"
# )
# 
# pair_sav_mh <- pairwise_sma(
#   model_formula = "sav_mean ~ mass_to_height", 
#   trait1 = sav_mean,
#   trait2 = mass_to_height,
#   data = "log"
# )
# 
# # fix this later with new mass data entry, though not sure why things worked with 
# # other models...
# pair_thickness_mh <- pairwise_sma(
#   model_formula = "thickness_mm_mean ~ mass_to_height",
#   trait1 = thickness_mm_mean,
#   trait2 = mass_to_height,
#   data = "log"
# )
# 
# sma(mass_to_height ~ thickness_mm_mean,
#     data = log_ind_traits %>% 
#       filter(mass_to_height != "-Inf"))
# 
# ggplot(data = log_ind_traits %>% 
#          filter(mass_to_height != "-Inf"),
#        aes(x = mass_to_height,
#            y = thickness_mm_mean)) +
#   geom_point(alpha = 0.5,
#              shape = 21) +
#   stat_ma_line(method = "SMA") +
#   # labs(title = model_formula) +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   labs(x = "log mass:height ratio",
#        y = "log thickness") +
#   coord_cartesian(ylim = c(-2.5, 0.6)) +
#   annotate(geom = "text", x = -1.5, y = -2.5, label = "R\U00B2 = 0.20, p < 0.001")
# 
# pair_thickness_sav <- pairwise_sma(
#   model_formula = "thickness_mm_mean ~ sav_mean", 
#   trait1 = thickness_mm_mean,
#   trait2 = sav_mean,
#   data = "log"
# )
# 
# pair_tdmc_height <- pairwise_sma(
#   model_formula = "tdmc_mean ~ maximum_height", 
#   trait1 = tdmc_mean,
#   trait2 = maximum_height,
#   data = "log"
# )
# 
# pair_tdmc_mh <- pairwise_sma(
#   model_formula = "tdmc_mean ~ mass_to_height", 
#   trait1 = tdmc_mean,
#   trait2 = mass_to_height,
#   data = "log"
# )
# 
# pair_tdmc_thickness <- pairwise_sma(
#   model_formula = "tdmc_mean ~ thickness_mm_mean", 
#   trait1 = tdmc_mean,
#   trait2 = thickness_mm_mean,
#   data = "log"
# )
# 
# pair_sta_height <- pairwise_sma(
#   model_formula = "sta_mean ~ maximum_height", 
#   trait1 = sta_mean,
#   trait2 = maximum_height,
#   data = "log"
# )
# 
# 
# 
# 
# pair_sta_thickness <- pairwise_sma(
#   model_formula = "sta_mean ~ thickness_mm_mean", 
#   trait1 = sta_mean,
#   trait2 = thickness_mm_mean,
#   data = "log"
# )
# 
# pair_sta_tdmc <- pairwise_sma(
#   model_formula = "sta_mean ~ tdmc_mean", 
#   trait1 = sta_mean,
#   trait2 = tdmc_mean,
#   data = "log"
# )
# 
# pair_sap_height <- pairwise_sma(
#   model_formula = "sap_mean ~ maximum_height", 
#   trait1 = sap_mean,
#   trait2 = maximum_height,
#   data = "log"
# )
# 
# pair_sap_sav <- pairwise_sma(
#   model_formula = "sap_mean ~ sav_mean", 
#   trait1 = sap_mean,
#   trait2 = sav_mean,
#   data = "log"
# )
# 
# pair_sap_thickness <- pairwise_sma(
#   model_formula = "sap_mean ~ thickness_mm_mean", 
#   trait1 = sap_mean,
#   trait2 = thickness_mm_mean,
#   data = "log"
# )
# 
# pair_sap_tdmc <- pairwise_sma(
#   model_formula = "sap_mean ~ tdmc_mean", 
#   trait1 = sap_mean,
#   trait2 = tdmc_mean,
#   data = "log"
# )
# 
# pair_sap_sta <- pairwise_sma(
#   model_formula = "sap_mean ~ sta_mean", 
#   trait1 = sap_mean,
#   trait2 = sta_mean,
#   data = "log"
# )     
# 
# sma(mass_to_height ~ sav_mean,
#         data = log_ind_traits)
# 
# ggplot(data = log_ind_traits,
#        aes(x = mass_to_height,
#            y = sav_mean)) +
#   geom_point(alpha = 0.5,
#              shape = 21) +
#   stat_ma_line(method = "SMA") +
#   # labs(title = model_formula) +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   labs(x = "log mass:height ratio",
#        y = "log surface area:volume ratio") +
#   coord_cartesian(ylim = c(4, 9.5)) +
#   annotate(geom = "text", x = -1.5, y = 4, label = "R\U00B2 = 0.049, p < 0.001")
# 
# pair_fvfm_height <- pairwise_sma(
#   model_formula = "fvfm_mean ~ maximum_height", 
#   trait1 = fvfm_mean,
#   trait2 = maximum_height,
#   data = "log"
# )     
# 
# pair_fvfm_sav <- pairwise_sma(
#   model_formula = "fvfm_mean ~ sav_mean", 
#   trait1 = fvfm_mean,
#   trait2 = sav_mean,
#   data = "log"
# )    
# 
# pair_fvfm_tdmc <- pairwise_sma(
#   model_formula = "fvfm_mean ~ tdmc_mean", 
#   trait1 = fvfm_mean,
#   trait2 = tdmc_mean,
#   data = "log"
# )
# 
# pair_fvfm_sta <- pairwise_sma(
#   model_formula = "fvfm_mean ~ sta_mean", 
#   trait1 = fvfm_mean,
#   trait2 = sta_mean,
#   data = "log"
# )
# 
# pair_fvfm_sap <- pairwise_sma(
#   model_formula = "fvfm_mean ~ sap_mean", 
#   trait1 = fvfm_mean,
#   trait2 = sap_mean,
#   data = "log"
# )
# 
# pair_ar_height <- pairwise_sma(
#   model_formula = "aspect_ratio_mean ~ maximum_height", 
#   trait1 = aspect_ratio_mean,
#   trait2 = maximum_height,
#   data = "log"
# )
# 
# # this also doesn't work!
# # pair_ar_mh <- pairwise_sma(
# #   model_formula = "aspect_ratio_mean ~ mass_to_height",
# #   trait1 = aspect_ratio_mean,
# #   trait2 = mass_to_height,
# #   data = "log"
# # )
# 
# pair_ar_sav <- pairwise_sma(
#   model_formula = "aspect_ratio_mean ~ sav_mean", 
#   trait1 = aspect_ratio_mean,
#   trait2 = sav_mean,
#   data = "log"
# )
# 
# pair_ar_thickness <- pairwise_sma(
#   model_formula = "aspect_ratio_mean ~ thickness_mm_mean", 
#   trait1 = aspect_ratio_mean,
#   trait2 = thickness_mm_mean,
#   data = "log"
# )
# 
# pair_ar_tdmc <- pairwise_sma(
#   model_formula = "aspect_ratio_mean ~ tdmc_mean", 
#   trait1 = aspect_ratio_mean,
#   trait2 = tdmc_mean,
#   data = "log"
# )
# 
# pair_ar_fvfm <- pairwise_sma(
#   model_formula = "aspect_ratio_mean ~ fvfm_mean", 
#   trait1 = aspect_ratio_mean,
#   trait2 = fvfm_mean,
#   data = "log"
# )
# 
# ar_plots <- (pair_ar_height[[3]] | pair_ar_fvfm[[3]] | pair_ar_sav[[3]] | pair_ar_tdmc[[3]] | pair_ar_thickness[[3]]) 
# 
# fvfm_plots <- 
#   (pair_fvfm_height[[3]] | pair_fvfm_sap[[3]] | pair_fvfm_sav[[3]] | pair_fvfm_tdmc[[3]] | pair_fvfm_sta[[3]] ) 
# 
# ar_plots/fvfm_plots
# 
# mh_row_scaled <- plot_grid(
#   NULL,
#   NULL,
#   NULL,
#   NULL,
#   NULL,
#   NULL,
#   NULL,
#   NULL,
#   nrow = 1
# )
# 
# mh_row <- plot_grid(
#   NULL,
#   NULL,
#   NULL,
#   NULL,
#   NULL,
#   NULL,
#   NULL,
#   NULL,
#   nrow = 1
# )
# 
# sav_row_scaled <- plot_grid(
#   NULL,
#   pair_sav_mh[[3]],
#   NULL,
#   NULL,
#   NULL,
#   NULL,
#   NULL,
#   NULL,
#   nrow = 1
# )
# 
# sav_row <- plot_grid(
#   pair_sav_height[[3]],
#   pair_sav_mh[[3]],
#   NULL,
#   NULL,
#   NULL,
#   NULL,
#   NULL,
#   NULL,
#   nrow = 1
# )
# 
# 
# 
# thickness_row_scaled <- plot_grid(
#   NULL,
#   pair_thickness_mh[[3]],
#   pair_thickness_sav[[3]],
#   NULL,
#   NULL,
#   NULL,
#   NULL,
#   NULL,
#   nrow = 1
# )
# 
# thickness_row <- plot_grid(
#   NULL,
#   NULL,
#   pair_thickness_sav[[3]],
#   NULL,
#   NULL,
#   NULL,
#   NULL,
#   NULL,
#   nrow = 1
# )
# 
# tdmc_row_scaled <- plot_grid(
#   NULL,
#   pair_tdmc_mh[[3]],
#   NULL,
#   pair_tdmc_thickness[[3]],
#   NULL,
#   NULL,
#   NULL,
#   NULL,
#   nrow = 1
# )
# 
# tdmc_row <- plot_grid(
#   pair_tdmc_height[[3]],
#   pair_tdmc_mh[[3]],
#   NULL,
#   pair_tdmc_thickness[[3]],
#   NULL,
#   NULL,
#   NULL,
#   NULL,
#   nrow = 1
# )
# 
# sta_row_scaled <- plot_grid(
#   NULL,
#   pair_sta_mh[[3]],
#   pair_sta_sav[[3]],
#   pair_sta_thickness[[3]],
#   pair_sta_tdmc[[3]],
#   NULL,
#   NULL,
#   NULL,
#   nrow = 1
# )
# 
# sta_row <- plot_grid(
#   pair_sta_height[[3]],
#   pair_sta_mh[[3]],
#   pair_sta_sav[[3]],
#   pair_sta_thickness[[3]],
#   pair_sta_tdmc[[3]],
#   NULL,
#   NULL,
#   NULL,
#   nrow = 1
# )
# 
# sap_row_scaled <- plot_grid(
#   NULL,
#   NULL,
#   pair_sap_sav[[3]],
#   pair_sap_thickness[[3]],
#   pair_sap_tdmc[[3]],
#   pair_sap_sta[[3]],
#   NULL,
#   NULL,
#   nrow = 1
# )
# 
# sap_row <- plot_grid(
#   pair_sap_height[[3]],
#   NULL,
#   pair_sap_sav[[3]],
#   NULL,
#   pair_sap_tdmc[[3]],
#   pair_sap_sta[[3]],
#   NULL,
#   NULL,
#   nrow = 1
# )
# 
# fvfm_row_scaled <- plot_grid(
#   pair_fvfm_height[[3]],
#   NULL,
#   pair_fvfm_sav[[3]],
#   NULL,
#   pair_fvfm_tdmc[[3]],
#   pair_fvfm_sta[[3]],
#   pair_fvfm_sap[[3]],
#   NULL,
#   nrow = 1
# )
# 
# fvfm_row <- plot_grid(
#   pair_fvfm_height[[3]],
#   NULL,
#   pair_fvfm_sav[[3]],
#   NULL,
#   pair_fvfm_tdmc[[3]],
#   pair_fvfm_sta[[3]],
#   pair_fvfm_sap[[3]],
#   NULL,
#   nrow = 1
# )
# 
# ar_row_scaled <- plot_grid(
#   NULL,
#   NULL,
#   pair_ar_sav[[3]], 
#   pair_ar_thickness[[3]],
#   NULL,
#   NULL,
#   NULL,
#   NULL,
#   nrow = 1
# )
# 
# ar_row <- plot_grid(
#   pair_ar_height[[3]],
#   NULL,
#   pair_ar_sav[[3]], 
#   pair_ar_thickness[[3]],
#   pair_ar_tdmc[[3]],
#   NULL,
#   NULL,
#   pair_ar_fvfm[[3]],
#   nrow = 1
# )
# 
# sma_grid_scaled <- plot_grid(mh_row_scaled,
#                              sav_row_scaled,
#                              thickness_row_scaled,
#                              tdmc_row_scaled,
#                              sta_row_scaled,
#                              sap_row_scaled,
#                              fvfm_row_scaled,
#                              ar_row_scaled,
#                              nrow = 8)
# 
# sma_grid <- plot_grid(mh_row,
#                       sav_row,
#                       thickness_row,
#                       tdmc_row,
#                       sta_row,
#                       sap_row,
#                       fvfm_row,
#                       ar_row,
#                       nrow = 8)
# 
# cowplot::save_plot(here::here(
#   "figures",
#   "tradeoffs",
#   paste0("sma-grid_", today(), ".jpg")),
#   base_height = 10,
#   sma_grid_scaled)
# 
# ggsave(filename = here::here(
#   "figures",
#   "tradeoffs",
#   paste0("sma-grid_", today(), ".jpg")),
#   plot = sma_grid,
#   width = 18,
#   height = 18,
#   units = "cm",
#   dpi = 300
# )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# -------------------- 5. Principal Components Analysis -------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section includes a Principal Components Analysis (PCA) to visualize the
# multivariate relationships between traits.

# I got a lot of help from David Zeleny's blog post: 
# https://www.davidzeleny.net/anadat-r/doku.php/en:pca_examples).


# ⟞ a. full model ---------------------------------------------------------

# ⟞ ⟞ i. PCA --------------------------------------------------------------

# trait by species PCA
pca_full <- rda(pca_mat_log, scale = TRUE)

# create a screeplot to visualize axis contributions
screeplot(pca_full, bstick = TRUE)

# look at the summary
summary(pca_full)
prcomp(pca_mat_log, center = TRUE, scale. = TRUE) %>% summary()

# proportion variance explained for downstream figure making
prop_PC1_full <- "62.6%"
prop_PC2_full <- "21.9%"

sp_permanova <- adonis2(pca_mat_log ~ sp_code, 
                        data = ind_traits_filtered,
                        method = "euclidean")
sp_permanova

adonis_pairwise <- pairwise.adonis2(pca_mat_log ~ sp_code, 
                 data = ind_traits_filtered)
rvam_pairwise <- pairwise.perm.manova(resp = pca_mat_log,
                     fact = ind_traits_filtered$sp_code)

# ⟞ ⟞ ii. loadings --------------------------------------------------------

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

# ⟞ ⟞ iii. PCA plots ------------------------------------------------------

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
                                  levels = algae_proposal_sciname_factors)) %>% 
  left_join(., fvfm_ind, by = "specimen_ID")

# trait vectors
PCAvect_full <- scores(pca_full, 
                       display = "species", 
                       choices = c(1, 2)) %>% 
  as.data.frame()

# plot PCA
plot_PCA_12_full <- ggplot() +
  PCA_aesthetics +
  geom_point(data = PCAscores_full, 
             aes(x = PC1, 
                 y = PC2#, 
                 # color = scientific_name, 
                 # shape = scientific_name,
                 # size = fvfm_mean
                 ) ,
             alpha = 0.7,
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
  scale_x_continuous(limits = c(-2, 3)) +
  scale_y_continuous(limits = c(-2.75, 2.25)) +
  PCA_theme() +
  labs(x = paste0("PC1 (", prop_PC1_full, ")"),
       y = paste0("PC2 (", prop_PC2_full, ")"),
       title = "(a) PC1 and PC2",
       subtitle = "BO, CC, CO, BF, DP, LAFA, PTCA, R, CYOS", 
       color = "Scientific name") 
plot_PCA_12_full

ggsave(here::here("figures",
                  "ordination",
                  paste("PCA-log_scale_full-model_", today(), ".jpg", sep = "")),
       plot_PCA_12_full,
       width = 12, height = 12, units = "cm", dpi = 300)

# ⟞ ⟞ iv. trait contributions to axes -------------------------------------

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
  mutate(contrib = (values/component_total)*100)

# The red dashed line on the graph above indicates the expected average contribution. If the contribution of the variables were uniform, the expected value would be 1/length(variables) = 1/10 = 10%
expected_average_full <- (1/length(unique(varcoord_full$trait)))*100

contrib_aesthetics_full <- list(
  geom_col(color = "black"),
  geom_vline(xintercept = expected_average_full,
             color = "#96a39e",
             linetype = 2),
  scale_fill_manual(values = trait_color_palette),
  scale_x_continuous(expand = expansion(c(0, 0.05)),
                     position = "top")
)

pc1_contrib_full <- varcoord_full %>% 
  filter(axis == "PC1") %>% 
  ggplot(aes(y = reorder(trait, contrib),
             x = contrib, 
             fill = trait)) +
  contrib_aesthetics_full +
  contrib_theme() +
  labs(title = "(b) Trait % contributions to PC1")

pc1_contrib_full

pc2_contrib_full <- varcoord_full %>% 
  filter(axis == "PC2") %>% 
  ggplot(aes(y = reorder(trait, contrib),
             x = contrib, 
             fill = trait)) +
  contrib_aesthetics_full +
  contrib_theme() +
  labs(title = "(c) Trait % contributions to PC2")

pc2_contrib_full

# scaled layout
contrib_together_full <- (plot_PCA_12_full) | ((pc1_contrib_full / pc2_contrib_full) + plot_layout(axis_titles = "collect")) 

contrib_together_full <- (free(plot_PCA_12_full) / 
                            (free(pc1_contrib_full) | free(pc2_contrib_full))) + 
  plot_layout(heights = c(1, 0.8))

# contrib_together_full <- ((free(plot_PCA_12_full) / (free(pc1_contrib_full)) | (free(pc2_contrib_full) / plot_spacer()))) + 
#   plot_layout(heights = c(1, 0.5, 0.5))

ggsave(here::here(
  "figures",
  "ordination",
  paste0("contributions_scale_full-model_", today(), ".jpg")),
  contrib_together_full,
  width = 14,
  height = 22,
  units = "cm",
  dpi = 300
)

# ⟞ b. reduced model ------------------------------------------------------

# ⟞ ⟞ i. PCA --------------------------------------------------------------

reduced_mat <- pca_mat_log %>% 
  select(`Dry:wet weight`, `Height:wet weight`, `Height`, `Thickness`)

# trait by species PCA
pca_reduced <- rda(reduced_mat, scale = TRUE)

# create a screeplot to visualize axis contributions
screeplot(pca_reduced, bstick = TRUE)

# look at the summary
summary(pca_reduced)

# proportion variance explained for downstream figure making
prop_PC1_reduced <- "48.1%"
prop_PC2_reduced <- "34.8%"

adonis_pairwise_reduced <- pairwise.adonis2(reduced_mat ~ sp_code, 
                                    data = ind_traits_filtered)
rvam_pairwise_reduced <- pairwise.perm.manova(reduced_mat,
                                              fact = ind_traits_filtered$sp_code)
rvam_pairwise_reduced 

# ⟞ ⟞ ii. loadings --------------------------------------------------------

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

# ⟞ ⟞ iii. PCA plots ------------------------------------------------------

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

# ⟞ ⟞ iv. trait contributions to axes -------------------------------------

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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---------------------- 6. species collection table ----------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

total_sample_collection_table <- ind_traits_filtered %>% 
  # rownames_to_column("specimen_ID") %>% 
  # select(specimen_ID) %>% 
  # left_join(., metadata_ind, by = "specimen_ID") %>% 
  # left_join(., coarse_traits, by = "sp_code") %>% 
  # filter(sp_code %in% algae_proposal) %>% 
  select(specimen_ID, scientific_name, site) %>% 
  group_by(scientific_name, site) %>% 
  count() %>% 
  pivot_wider(names_from = "site",
              values_from = "n") %>% 
  select(scientific_name, Bullito, `Arroyo Quemado`, Naples, `Isla Vista`,
         Mohawk, Carpinteria) %>% 
  adorn_totals(c("row", "col")) %>% 
  rename(`Scientific name` = scientific_name) %>% 
  flextable() %>% 
  colformat_num(
    part = "all",
    na_str = "-"
  ) %>% 
  autofit() %>% 
  fit_to_width(7) %>% 
  font(fontname = "Times New Roman",
       part = "all")

total_sample_collection_table

total_sample_collection_table %>%
  save_as_docx(path = here::here(
    "tables",
    "sample-tables",
    paste0("total-samples_all-data_", today(), ".docx")
    ))

# total_samples_with_all_data <- ind_traits %>% 
#   filter(sp_code %in% algae_proposal) %>% 
#   select(scientific_name, site,
#          maximum_height, thickness_mm_mean, sta_mean, frond_dmc_mean,
#          sav_mean, sap_mean, aspect_ratio_mean, frond_length_mean,
#          frond_width_mean, fvfm_mean, total_wet, total_dry, total_volume,
#          total_dmc) %>% 
#   drop_na(maximum_height, thickness_mm_mean, sta_mean, frond_dmc_mean,
#           sav_mean, sap_mean, aspect_ratio_mean, frond_length_mean,
#           frond_width_mean, fvfm_mean, total_wet, total_dry, total_volume,
#           total_dmc) %>% 
#   select(scientific_name, site) %>% 
#   group_by(scientific_name, site) %>% 
#   count() %>% 
#   pivot_wider(names_from = "site",
#               values_from = "n") %>% 
#   select(scientific_name, Bullito, `Arroyo Quemado`, Naples, `Isla Vista`,
#          Mohawk, Carpinteria) %>% 
#   adorn_totals(c("row", "col")) %>% 
#   rename(`Scientific name` = scientific_name) %>% 
#   flextable() %>% 
#   colformat_num(
#     part = "all",
#     na_str = "-"
#   ) %>% 
#   autofit() %>% 
#   fit_to_width(10) %>% 
#   font(fontname = "Times New Roman",
#        part = "all")
# 
# total_samples_with_all_data 
# 
# total_samples_with_all_data %>%
#   save_as_docx(path = here::here(
#     "tables",
#     "sample-tables",
#     paste0("total-samples_all-data_", today(), ".docx")
#     ))
