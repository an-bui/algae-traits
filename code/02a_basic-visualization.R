
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

boxplots <- ind_traits_filtered %>% 
  pivot_longer(cols = maximum_height:sap_mean,
               names_to = "trait",
               values_to = "value") %>% 
  nest(.by = trait,
       data = everything()) %>% 
  mutate(units = map(
    trait,
    ~ case_match(
      .x,
      "maximum_height" ~ "Maximum height (cm)", 
      "thickness_mm_mean" ~ "Thickness (mm)", 
      "frond_area_scaled" ~ "Surface area (mm\U00B2)",
      "height_ww" ~ "Height:Wet weight (cm/wet mg)",
      "total_dmc" ~ "Dry:wet weight", 
      "height_vol" ~ "Height:volume (cm/mL)",
      "sav_scaled" ~ "Surface area:volume (mm\U00B2/mL)", 
      "sta_scaled" ~ "Surface area:dry weight (mm\U00B2/dry mg)",
      "sap_mean" ~ "Surface area:perimeter (mm\U00B2/mm)"
    )
  )) %>% 
  mutate(boxplot = map2(
    data, units,
    ~ .x %>% 
      mutate( 
      sp_code = fct_relevel(sp_code,
                                  "R", "BO", "CO", "CC", "BF", "DP", "LAFA", "PTCA", "CYOS")) %>% 
      ggplot(aes(x = sp_code,
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

# maximum height
pluck(boxplots, 4, 1)

# thickness
pluck(boxplots, 4, 2)

# surface area
pluck(boxplots, 4, 3)

# height:wet weight
pluck(boxplots, 4, 4)

# dry:wet weight
pluck(boxplots, 4, 5)

# height:volume
pluck(boxplots, 4, 6)

# SA:V
pluck(boxplots, 4, 7)

# surface area:dry weight
pluck(boxplots, 4, 8)

# SA:P
pluck(boxplots, 4, 9)


# ⟞ b. multipanel plot ----------------------------------------------------

boxplot_multipanel <- 
  (pluck(boxplots, 4, 1) + pluck(boxplots, 4, 2) + pluck(boxplots, 4, 3)) /
  (pluck(boxplots, 4, 4) + pluck(boxplots, 4, 5) + pluck(boxplots, 4, 6)) /
  (pluck(boxplots, 4, 7) + pluck(boxplots, 4, 8) + pluck(boxplots, 4, 9)) 

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

distributions <- ind_traits_filtered %>% 
  pivot_longer(cols = maximum_height:sap_mean,
               names_to = "trait",
               values_to = "value") %>% 
  nest(.by = trait,
       data = everything()) %>% 
  mutate(units = map(
    trait,
    ~ case_match(
      .x,
      "maximum_height" ~ "Maximum height (cm)", 
      "thickness_mm_mean" ~ "Thickness (mm)", 
      "frond_area_scaled" ~ "Surface area (mm\U00B2)",
      "height_ww" ~ "Height:Wet weight (cm/wet mg)",
      "total_dmc" ~ "Dry:wet weight", 
      "height_vol" ~ "Height:volume (cm/mL)",
      "sav_scaled" ~ "Surface area:volume (mm\U00B2/mL)", 
      "sta_scaled" ~ "Surface area:dry weight (mm\U00B2/dry mg)",
      "sap_mean" ~ "Surface area:perimeter (mm\U00B2/mm)"
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
      "maximum_height" ~ max_height_col, 
      "thickness_mm_mean" ~ thickness_col, 
      "frond_area_scaled" ~ sta_col,
      "height_ww" ~ mass_to_height_col,
      "total_dmc" ~ total_dmc_col, 
      "height_vol" ~ mass_to_height_col,
      "sav_scaled" ~ sav_col, 
      "sta_scaled" ~ sta_col,
      "sap_mean" ~ sap_col
    )
  )) %>% 
  mutate(distribution_plot = pmap(
    list(w = data, x = bins, y = colors, z = units),
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
    list(w = data, x = bins, y = colors, z = units),
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
    list(w = data, x = bins, y = colors, z = units),
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
    list(x = data, y = units),
    function(x, y) ggplot(data = x,
                          aes(sample = value)) +
      geom_qq_line(color = "darkgrey") +
      geom_qq(shape = 21) +
      labs(title = paste0(y, " (no transformation)")) +
      distribution_theme
  )) %>% 
  mutate(qq_log = pmap(
    list(x = data, y = units),
    function(x, y) ggplot(data = x,
                          aes(sample = log(value))) +
      geom_qq_line(color = "darkgrey") +
      geom_qq(shape = 21) +
      labs(title = paste0(y, " (log transform)")) +
      distribution_theme
  )) %>% 
  mutate(qq_sqrt = pmap(
    list(x = data, y = units),
    function(x, y) ggplot(data = x,
                          aes(sample = sqrt(value))) +
      geom_qq_line(color = "darkgrey") +
      geom_qq(shape = 21) +
      labs(title = paste0(y, " (square root transform)")) +
      distribution_theme
  ))

# max height
pluck(distributions, 7, 1) # right skewed
pluck(distributions, 10, 1)
pluck(distributions, 8, 1) # log transform looks normal (ish)
pluck(distributions, 11, 1)
pluck(distributions, 9, 1)
pluck(distributions, 12, 1)

# thickness
pluck(distributions, 7, 2) # right skewed
pluck(distributions, 10, 2)
pluck(distributions, 8, 2) # log transform looks normal
pluck(distributions, 11, 2)
pluck(distributions, 9, 2)
pluck(distributions, 12, 2)

# surface area
pluck(distributions, 7, 3) # right skewed
pluck(distributions, 10, 3)
pluck(distributions, 8, 3) # log transform looks normal
pluck(distributions, 11, 3)
pluck(distributions, 9, 3)
pluck(distributions, 12, 3)


# height:wet weight
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

# height:volume
pluck(distributions, 7, 6) # right skewed
pluck(distributions, 10, 6)
pluck(distributions, 8, 6) # log transform is normal
pluck(distributions, 11, 6) 
pluck(distributions, 9, 6)
pluck(distributions, 12, 6)


# surface area:volume
pluck(distributions, 7, 7) # right skewed
pluck(distributions, 10, 7)
pluck(distributions, 8, 7) 
pluck(distributions, 11, 7) # log transform is normal
pluck(distributions, 9, 7)
pluck(distributions, 12, 7)


# surface area:dry weight
pluck(distributions, 7, 8) 
pluck(distributions, 10, 8)
pluck(distributions, 8, 8) # log transform is normal
pluck(distributions, 11, 8)
pluck(distributions, 9, 8)
pluck(distributions, 12, 8)


# surface area:perimeter
pluck(distributions, 7, 9) 
pluck(distributions, 10, 9)
pluck(distributions, 8, 9) 
pluck(distributions, 11, 9) # normal
pluck(distributions, 9, 9)
pluck(distributions, 12, 9)

# transform all traits


# ⟞ b. multipanel plot ----------------------------------------------------

distributions_multipanel <- 
  (pluck(distributions, 8, 1) + pluck(distributions, 8, 2) + pluck(distributions, 8, 3)) /
  (pluck(distributions, 8, 4) + pluck(distributions, 8, 5) + pluck(distributions, 8, 6)) /
  (pluck(distributions, 8, 7) + pluck(distributions, 8, 8) + pluck(distributions, 8, 9))

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
