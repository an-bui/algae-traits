
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 0. source -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# only need to run this once per session
source(here::here("code", "01a_trait-cleaning.R"))

# basic theme
boxplot_theme <- list(
    theme_bw(),
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 12),
          plot.title = element_text(size = 26))
)

distribution_theme <- list(
  theme_bw(),
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 22))
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------ 1. boxplots ------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section includes code to generate a bunch of boxplots showing trait
# distributions across species. Each point is an individual, and the boxplots
# depict the median, IQR, and 1.5*IQR. The outliers aren't shown.

# ⟞ a. generating plots ---------------------------------------------------

boxplots <- ind_traits %>% 
  filter(sp_code %in% c("BO", "CC", "BF", "DP", "LAFA", "PTCA", "R", "CYOS")) %>% 
  select(specimen_ID, scientific_name, sp_code, pigment_type,
         fvfm_mean,
         maximum_height, 
         mass_to_height, total_dry, 
         total_dmc, total_wet,
         total_volume,
         sav_scaled, frond_area_scaled,
         thickness_mm_mean, 
         sta_scaled,
         sap_mean, frond_peri_scaled,
         aspect_ratio_mean, frond_length_scaled, frond_width_scaled) %>% 
  pivot_longer(cols = fvfm_mean:frond_width_scaled,
               names_to = "trait",
               values_to = "value") %>% 
  nest(.by = trait,
       data = everything()) %>% 
  mutate(trait = map(
    trait,
    ~ case_match(
      .x, 
      "fvfm_mean" ~ "fvfm",
      "maximum_height" ~ "max_height", 
      "mass_to_height" ~ "mass_to_height", 
      "total_dry" ~ "total_dry",
      "total_dmc" ~ "total_dmc",
      "total_wet" ~ "total_wet",
      "total_volume" ~ "total_volume",
      "sav_scaled" ~ "sav", 
      "frond_area_scaled" ~ "area", 
      "thickness_mm_mean" ~ "thickness", 
      "sta_scaled" ~ "sta", 
      "sap_mean" ~ "sap", 
      "frond_peri_scaled" ~ "peri",
      "aspect_ratio_mean" ~ "aspect_ratio", 
      "frond_length_scaled" ~ "length", 
      "frond_width_scaled" ~ "width"
    )
  )) %>% 
  mutate(units = map(
    trait,
    ~ case_match(
      .x,
      "fvfm" ~ "Fv/Fm",
      "max_height" ~ "Maximum height (cm)",
      "mass_to_height" ~ "Dry weight:height (dry mg/cm)",
      "total_dry" ~ "Individual dry weight (mg)",
      "total_dmc" ~ "Individual dry:wet weight",
      "total_wet" ~ "Individual wet weight (mg)",
      "total_volume" ~ "Volume (mL)",
      "sav" ~ "Surface area:volume ratio (mm\U00B2/mL)",
      "area" ~ "Surface area (mm\U00B2)", 
      "thickness" ~ "Thickness (mm)",
      "sta" ~ "Surface area:dry weight ratio (mm\U00B2/dry mg)",
      "sap" ~ "Surface area:perimeter ratio (mm\U00B2/mm)",
      "peri" ~ "Perimeter (mm)",
      "aspect_ratio" ~ "Aspect ratio",
      "length" ~ "Length (cm)", 
      "width" ~ "Width (cm)"
    )
  )) %>% 
  mutate(boxplot = map2(
    data, units,
    ~ .x %>% 
      mutate( 
      sp_code_label = case_when(
        sp_code == "Nandersoniana" ~ "NA",
        TRUE ~ sp_code
      ),
      sp_code_label = fct_relevel(sp_code_label,
                                  "R", "BO", "CC", "BF", "DP", "LAFA", "PTCA", "CYOS")) %>% 
      ggplot(aes(x = sp_code_label,
                 y = value)) +
      geom_boxplot(aes(fill = sp_code),
                   outliers = FALSE) +
      geom_point(alpha = 0.2,
                 position = position_jitter(width = 0.1,
                                            seed = 666),
                 shape = 21) +
      scale_fill_manual(values = algae_spcode_colors) +
      # facet_grid(cols = vars(pigment_type), 
      #            scales = "free_x") +
      scale_x_discrete(labels = scales::label_wrap(10)) +
      labs(title = .y) +
      boxplot_theme +
      theme(axis.title = element_blank(),
            axis.title.x = element_blank())
  ))

# FvFm
pluck(boxplots, 4, 1)

# maximum height
pluck(boxplots, 4, 2)

# mass:height
pluck(boxplots, 4, 3)

# individual dry weight
pluck(boxplots, 4, 4)

# individual dry matter content
pluck(boxplots, 4, 5)

# individual wet weight
pluck(boxplots, 4, 6)

# individual volume
pluck(boxplots, 4, 7)

# SA:V
pluck(boxplots, 4, 8)

# surface area
pluck(boxplots, 4, 9)

# Thickness
pluck(boxplots, 4, 10)

# surface area:dry mass
pluck(boxplots, 4, 11)

# sa:p
pluck(boxplots, 4, 12)

# perimeter
pluck(boxplots, 4, 13)

# aspect ratio
pluck(boxplots, 4, 14)

# frond length
pluck(boxplots, 4, 15)

# frond width
pluck(boxplots, 4, 16)


# ⟞ b. multipanel plot ----------------------------------------------------

boxplot_multipanel <- 
  # sa:dw, sa:V, SA
  (pluck(boxplots, 4, 11) + pluck(boxplots, 4, 8) + pluck(boxplots, 4, 9)) /
  # v, sa:p, p
  (pluck(boxplots, 4, 7) + pluck(boxplots, 4, 12) + pluck(boxplots, 4, 13)) /
  # aspect, l, w
  # (pluck(boxplots, 4, 14) + pluck(boxplots, 4, 15) + pluck(boxplots, 4, 16)) /
  # dw:ww, dw, ww
  (pluck(boxplots, 4, 5) + pluck(boxplots, 4, 4) + pluck(boxplots, 4, 6)) /
  # thickness, dw:h, max height
  (pluck(boxplots, 4, 10) + pluck(boxplots, 4, 3) + pluck(boxplots, 4, 2))

# ⟞ c. saving outputs -----------------------------------------------------

# traits in the boxplot data frame
boxplot_traits <- c("fvfm", "thickness", "volume", "sap", "sav", 
                    "max_height", "sta", "mass_to_height", "aspect_ratio", 
                    "chlA", "frond_dmc", "total_dmc")

# function to save boxplot
save_boxplot <- function(trait, boxplot) {
  ggsave(here("figures",
              "basic-visualizations",
              "boxplots",
              paste0("boxplot-", trait, "_", today(), ".jpg")),
         boxplot,
         width = 16, height = 10, units = "cm",
         dpi = 300)
}

# for loop to save boxplots
for(i in 1:length(boxplot_traits)) {
  
  # save the trait as an object
  trait <- boxplot_traits[i]
  # save the boxplot as an object
  boxplot <- pluck(boxplots, 4, i)
  
  # run the function to save the boxplot
  save_boxplot(trait = trait, 
               boxplot = boxplot)
  
}

# ggsave(here("figures",
#             "basic-visualizations",
#             "boxplots",
#             paste0("multipanel_boxplots_", today(), ".jpg")),
#        boxplot_multipanel,
#        width = 24,
#        height = 18,
#        units = "cm",
#        dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---------------------------- 2. distributions ---------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section includes code to generate histograms and qq plots for each trait
# in three forms: 1) untransformed, 2) log transformed, and 3) square root
# transformed. Because many of the traits followed a skewed distribution, we
# wanted to know which traits would be candidates for transformation to follow
# a normal distribution.

# ⟞ a. generating plots ---------------------------------------------------

distributions <- ind_traits %>% 
  filter(sp_code %in% c("BO", "CC", "BF", "DP", "LAFA", "PTCA", "R", "CYOS")) %>% 
  select(specimen_ID, scientific_name, sp_code, pigment_type,
         fvfm_mean,
         maximum_height, 
         mass_to_height, total_dry, 
         total_dmc, total_wet,
         total_volume,
         sav_scaled, frond_area_scaled,
         thickness_mm_mean, 
         sta_scaled,
         sap_mean, frond_peri_scaled,
         aspect_ratio_mean, frond_length_scaled, frond_width_scaled) %>% 
  pivot_longer(fvfm_mean:frond_width_scaled,
               names_to = "trait",
               values_to = "value") %>% 
  nest(.by = "trait",
       data = everything()) %>% 
  # shortening trait names
  mutate(trait = map(
    trait,
    ~ case_match(
      .x,
      "fvfm_mean" ~ "fvfm",
      "maximum_height" ~ "max_height", 
      "mass_to_height" ~ "mass_to_height", 
      "total_dry" ~ "total_dry",
      "total_dmc" ~ "total_dmc",
      "total_wet" ~ "total_wet",
      "total_volume" ~ "total_volume",
      "sav_scaled" ~ "sav", 
      "frond_area_scaled" ~ "area", 
      "thickness_mm_mean" ~ "thickness", 
      "sta_scaled" ~ "sta", 
      "sap_mean" ~ "sap", 
      "frond_peri_scaled" ~ "peri",
      "aspect_ratio_mean" ~ "aspect_ratio", 
      "frond_length_scaled" ~ "length", 
      "frond_width_scaled" ~ "width"
    )
  )) %>% 
  # full names of traits
  mutate(trait_name = map(
    trait,
    ~ case_match(
      .x,
      "fvfm" ~ "Fv/Fm",
      "max_height" ~ "Maximum height (cm)",
      "mass_to_height" ~ "Dry weight:height (dry mg/cm)",
      "total_dry" ~ "Individual dry weight (mg)",
      "total_dmc" ~ "Individual dry:wet weight",
      "total_wet" ~ "Individual wet weight (mg)",
      "total_volume" ~ "Volume (mL)",
      "sav" ~ "Surface area:volume ratio (mm\U00B2/mL)",
      "area" ~ "Surface area (mm\U00B2)", 
      "thickness" ~ "Thickness (mm)",
      "sta" ~ "Surface area:dry weight ratio (mm\U00B2/dry mg)",
      "sap" ~ "Surface area:perimeter ratio (mm\U00B2/mm)",
      "peri" ~ "Perimeter (mm)",
      "aspect_ratio" ~ "Aspect ratio",
      "length" ~ "Length (cm)", 
      "width" ~ "Width (cm)"
    )
  )) %>% 
  mutate(length = map(
    data,
    ~ .x %>% 
      drop_na(value) %>% 
      nrow(.)
  )) %>% 
  mutate(bins = map(
    length,
    ~ round((.x^(1/3))*2)
  )) %>% 
  # colors for figures
  mutate(colors = map(
    trait,
    ~ case_match(
      .x,
      # 
      # sta_col <- "#BD973D"
      # sav_col <- "#CC7556"
      # sap_col <- "#CC7556"
      # aspect_ratio_col <- "#CC7556"
      # frond_dmc_col <- "#BE5A47"
      # total_dmc_col <- "#BE5A47"
      # thickness_col <- "#BE5A47"
      # mass_to_height_col <- "#6B6D9F"
      # max_height_col <- "#4D5B75"
      # fvfm_col <- "#b2d8d8"
      # volume_col <- "#BD973D"
      "fvfm" ~ fvfm_col,
      "max_height" ~ max_height_col,
      "mass_to_height" ~ mass_to_height_col,
      "total_dry" ~ total_dmc_col,
      "total_dmc" ~ total_dmc_col,
      "total_wet" ~ total_dmc_col,
      "total_volume" ~ volume_col,
      "sav" ~ sav_col,
      "area" ~ sav_col, 
      "thickness" ~ thickness_col,
      "sta" ~ sta_col,
      "sap" ~ sap_col,
      "peri" ~ sap_col,
      "aspect_ratio" ~ aspect_ratio_col,
      "length" ~ aspect_ratio_col, 
      "width" ~ aspect_ratio_col
    )
  )) %>% 
  mutate(distribution_plot = pmap(
    list(w = data, x = bins, y = colors, z = trait_name),
    function(w, x, y, z) ggplot(data = w,
                             aes(x = value)) +
      geom_histogram(bins = x,
                     fill = y,
                     color = "black") +
      labs(title = paste0(z, " (no transformation)"),
           sample = "Value",
           y = "Count") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      distribution_theme
  )) %>% 
  mutate(distribution_plot_log = pmap(
    list(w = data, x = bins, y = colors, z = trait_name),
    function(w, x, y, z) ggplot(data = w,
                             aes(x = log(value))) +
      geom_histogram(bins = x,
                     fill = y,
                     color = "black") +
      labs(title = paste0(z, " (log transform)"),
           sample = "natural log(value)",
           y = "Count") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      distribution_theme
  )) %>% 
  mutate(distribution_plot_sqrt = pmap(
    list(w = data, x = bins, y = colors, z = trait_name),
    function(w, x, y, z) ggplot(data = w,
                             aes(x = sqrt(value))) +
      geom_histogram(bins = x,
                     fill = y,
                     color = "black") +
      labs(title = paste0(z, " (square root transform)"),
           sample = "square root (value)",
           y = "Count") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      distribution_theme
  )) %>% 
  mutate(qq = pmap(
    list(x = data, y = trait),
    function(x, y) ggplot(data = x,
                          aes(sample = value)) +
      geom_qq_line(color = "darkgrey") +
      geom_qq(shape = 21) +
      labs(title = paste0(y, " (no transformation)")) +
      distribution_theme
  )) %>% 
  mutate(qq_log = pmap(
    list(x = data, y = trait),
    function(x, y) ggplot(data = x,
                          aes(sample = log(value))) +
      geom_qq_line(color = "darkgrey") +
      geom_qq(shape = 21) +
      labs(title = paste0(y, " (log transform)")) +
      distribution_theme
  )) %>% 
  mutate(qq_sqrt = pmap(
    list(x = data, y = trait),
    function(x, y) ggplot(data = x,
                          aes(sample = sqrt(value))) +
      geom_qq_line(color = "darkgrey") +
      geom_qq(shape = 21) +
      labs(title = paste0(y, " (square root transform)")) +
      distribution_theme
  ))

# fvfm
pluck(distributions, 7, 1) # right skewed
pluck(distributions, 10, 1)
pluck(distributions, 8, 1) # log transform looks normal (ish)
pluck(distributions, 11, 1)
pluck(distributions, 9, 1)
pluck(distributions, 12, 1)

# max height
pluck(distributions, 7, 2) # right skewed
pluck(distributions, 10, 2)
pluck(distributions, 8, 2) # log transform looks normal
pluck(distributions, 11, 2)
pluck(distributions, 9, 2)
pluck(distributions, 12, 2)

# dry weight:height
pluck(distributions, 7, 3) # right skewed
pluck(distributions, 10, 3)
pluck(distributions, 8, 3) # log transform looks normal
pluck(distributions, 11, 3)
pluck(distributions, 9, 3)
pluck(distributions, 12, 3)

# total dry weight
pluck(distributions, 7, 4) # right skewed
pluck(distributions, 10, 4)
pluck(distributions, 8, 4) # log transform looks normal
pluck(distributions, 11, 4)
pluck(distributions, 9, 4)
pluck(distributions, 12, 4)

# dry:wet weight
pluck(distributions, 7, 5) # right skewed
pluck(distributions, 10, 5)
pluck(distributions, 8, 5) # log transform looks better
pluck(distributions, 11, 5)
pluck(distributions, 9, 5)
pluck(distributions, 12, 5)

# total wet weight
pluck(distributions, 7, 6) # right skewed
pluck(distributions, 10, 6)
pluck(distributions, 8, 6) # log transform is normal
pluck(distributions, 11, 6) 
pluck(distributions, 9, 6)
pluck(distributions, 12, 6)

# individual volume
pluck(distributions, 7, 7) # right skewed
pluck(distributions, 10, 7)
pluck(distributions, 8, 7) 
pluck(distributions, 11, 7) # log transform is normal
pluck(distributions, 9, 7)
pluck(distributions, 12, 7)

# surface area:volume
pluck(distributions, 7, 8) 
pluck(distributions, 10, 8)
pluck(distributions, 8, 8) # log transform is normal
pluck(distributions, 11, 8)
pluck(distributions, 9, 8)
pluck(distributions, 12, 8)
# no transformations look good

# surface area
pluck(distributions, 7, 9) 
pluck(distributions, 10, 9)
pluck(distributions, 8, 9) 
pluck(distributions, 11, 9) # normal
pluck(distributions, 9, 9)
pluck(distributions, 12, 9)

# thickness
pluck(distributions, 7, 10) # right skewed
pluck(distributions, 10, 10)
pluck(distributions, 8, 10) 
pluck(distributions, 11, 10) # log transformation looks better
pluck(distributions, 9, 10)
pluck(distributions, 12, 10)

# frond DMC
pluck(distributions, 7, 11) # right skewed
pluck(distributions, 10, 11)
pluck(distributions, 8, 11) # log transform looks normal
pluck(distributions, 11, 11)
pluck(distributions, 9, 11)
pluck(distributions, 12, 11)
# no transformations look good

# surface area:perimeter
pluck(distributions, 7, 12) # right skewed
pluck(distributions, 10, 12)
pluck(distributions, 8, 12) # log transform looks normal
pluck(distributions, 11, 12)
pluck(distributions, 9, 12)
pluck(distributions, 12, 12)

# perimeter
pluck(distributions, 7, 13) # right skewed
pluck(distributions, 10, 13)
pluck(distributions, 8, 13) # log transform looks normalish
pluck(distributions, 11, 13)
pluck(distributions, 9, 13)
pluck(distributions, 12, 13)

# aspect ratio
pluck(distributions, 7, 14)
pluck(distributions, 10, 14)
pluck(distributions, 8, 14) # log transform looks normal
pluck(distributions, 11, 14)
pluck(distributions, 9, 14)
pluck(distributions, 12, 14)

# length
pluck(distributions, 7, 15) # right skewed
pluck(distributions, 10, 15)
pluck(distributions, 8, 15) # log transform looks normal
pluck(distributions, 11, 15)
pluck(distributions, 9, 15)
pluck(distributions, 12, 15)

# width
pluck(distributions, 7, 16) # right skewed
pluck(distributions, 10, 16)
pluck(distributions, 8, 16) # log transform looks normal
pluck(distributions, 11, 16)
pluck(distributions, 9, 16)
pluck(distributions, 12, 16)

# transform all traits


# ⟞ b. multipanel plot ----------------------------------------------------

distributions_multipanel <- 
  # sa:dw, sa:V, SA
  (pluck(distributions, 8, 11) + pluck(distributions, 8, 8) + pluck(distributions, 8, 9)) /
  # v, sa:p, p
  (pluck(distributions, 8, 7) + pluck(distributions, 8, 12) + pluck(distributions, 8, 13)) /
  # aspect, l, w
  # (pluck(distributions, 8, 14) + pluck(distributions, 8, 15) + pluck(distributions, 8, 16)) /
  # dw:ww, dw, ww
  (pluck(distributions, 8, 5) + pluck(distributions, 8, 4) + pluck(distributions, 8, 6)) /
  # thickness, dw:h, max height
  (pluck(distributions, 8, 10) + pluck(distributions, 8, 3) + pluck(distributions, 8, 2))

# ⟞ c. saving outputs -----------------------------------------------------

# traits in the distributions data frame

# "maximum_height" ~    "max_height",
# "mass_to_height" ~    "mass_to_height",
# "sav_mean" ~          "sav",
# "thickness_mm_mean" ~ "thickness",
# "tdmc_mean" ~         "TDMC",
# "sta_mean" ~          "sta",
# "sap_mean" ~          "sap",
# "fvfm_mean" ~         "fvfm",
# "aspect_ratio_mean" ~ "aspect_ratio",
# "frond_length_mean" ~ "thallus_length",
# "frond_width_mean" ~  "thallus_width",
# "total_wet" ~         "weight_wet",
# "total_dry" ~         "weight_dry",
# "total_volume" ~      "volume",
# "chlA_mean" ~         "chlA"
distributions_traits <- c("max_height",
                          "mass_to_height",
                          "sav",
                          "thickness",
                          "frond_dmc",
                          "sta",
                          "sap",
                          "fvfm",
                          "aspect_ratio",
                          "thallus_length",
                          "thallus_width",
                          "weight_wet",
                          "weight_dry",
                          "volume",
                          "chlA",
                          "total_dmc")

# function to save histograms
save_hist_qq <- function(trait, 
                         no_transform, 
                         log_transform, 
                         sqrt_transform,
                         no_transform_qq,
                         log_transform_qq,
                         sqrt_transform_qq) {
  ggsave(here("figures",
              "basic-visualizations",
              "distributions",
              paste0(trait, "_hist_no-transform_", today(), ".jpg")),
         no_transform,
         width = 12, height = 8, units = "cm",
         dpi = 200)
  
  ggsave(here("figures",
              "basic-visualizations",
              "distributions",
              paste0(trait, "_hist_log-transform_", today(), ".jpg")),
         log_transform,
         width = 12, height = 8, units = "cm",
         dpi = 200)
  
  ggsave(here("figures",
              "basic-visualizations",
              "distributions",
              paste0(trait, "_hist_sqrt-transform_", today(), ".jpg")),
         sqrt_transform,
         width = 12, height = 8, units = "cm",
         dpi = 200)
  
  ggsave(here("figures",
              "basic-visualizations",
              "distributions",
              paste0(trait, "_qq_no-transform_", today(), ".jpg")),
         no_transform_qq,
         width = 12, height = 8, units = "cm",
         dpi = 200)
  
  ggsave(here("figures",
              "basic-visualizations",
              "distributions",
              paste0(trait, "_qq_log-transform_", today(), ".jpg")),
         log_transform_qq,
         width = 12, height = 8, units = "cm",
         dpi = 200)
  
  ggsave(here("figures",
              "basic-visualizations",
              "distributions",
              paste0(trait, "_qq_sqrt-transform_", today(), ".jpg")),
         sqrt_transform_qq,
         width = 12, height = 8, units = "cm",
         dpi = 200)
}

# for loop to save boxplots
for(i in 1:length(distributions_traits)) {
  
  # save the trait as an object
  trait <- distributions_traits[i]
  # save the histograms as objects
  no_transform <- pluck(distributions, 5, i)
  log_transform <- pluck(distributions, 6, i)
  sqrt_transform <- pluck(distributions, 7, i)
  # save the qq plots as objects
  no_transform_qq <- pluck(distributions, 8, i)
  log_transform_qq <- pluck(distributions, 9, i)
  sqrt_transform_qq <- pluck(distributions, 10, i)
  
  # run the function to save the histograms
  save_hist_qq(trait = trait,
               no_transform = no_transform,
               log_transform = log_transform,
               sqrt_transform = sqrt_transform,
               no_transform_qq = no_transform_qq,
               log_transform_qq = log_transform_qq,
               sqrt_transform_qq = sqrt_transform_qq)
  
}

# ggsave(here("figures",
#             "basic-visualizations",
#             "distributions",
#             paste0("multipanel_hist_", today(), ".jpg")),
#        distributions_multipanel,
#        width = 24,
#        height = 18,
#        units = "cm",
#        dpi = 300)
