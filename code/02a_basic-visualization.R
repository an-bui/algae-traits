
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------- 0. source -------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# only need to run this once per session
source(here::here("code", "01a_trait-cleaning.R"))

# basic theme
basic_theme <- list(
    theme_bw(),
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.title = element_text(size = 14))
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------ 1. boxplots ------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This section includes code to generate a bunch of boxplots showing trait
# distributions across species. Each point is an individual, and the boxplots
# depict the median, IQR, and 1.5*IQR. The outliers aren't shown.

# ⟞ a. generating plots ---------------------------------------------------

boxplots <- ind_traits %>% 
  filter(sp_code %in% algae_proposal) %>% 
  select(specimen_ID, scientific_name, sp_code, pigment_type,
         fvfm_mean, thickness_mm_mean, total_volume,
         sap_mean, sav_mean, maximum_height, sta_mean, mass_to_height,
         aspect_ratio_mean, chlA_mean, frond_dmc_mean, total_dmc) %>% 
  pivot_longer(cols = fvfm_mean:total_dmc,
               names_to = "trait",
               values_to = "value") %>% 
  nest(.by = trait,
       data = everything()) %>% 
  mutate(trait = map(
    trait,
    ~ case_match(
      .x, 
      "fvfm_mean" ~ "fvfm",
      "thickness_mm_mean" ~ "thickness",
      "total_volume" ~ "volume",
      "sap_mean" ~ "sap",
      "sav_mean" ~ "sav",
      "maximum_height" ~ "max_height",
      "sta_mean" ~ "sta",
      "mass_to_height" ~ "mass_to_height",
      "aspect_ratio_mean" ~ "aspect_ratio",
      "chlA_mean" ~ "chlA",
      "frond_dmc_mean" ~ "frond_dmc",
      "total_dmc" ~ "total_dmc"
    )
  )) %>% 
  mutate(units = map(
    trait,
    ~ case_match(
      .x,
      "fvfm" ~ "Fv/Fm",
      "thickness" ~ "Mean thickness (mm)",
      "volume" ~ "Volume (mL)",
      "sap" ~ "Surface area:perimeter (mm\U00B2/mm)",
      "sav" ~ "Surface area:volume (mm\U00B2/mL)",
      "max_height" ~ "Maximum height (cm)",
      "sta" ~ "Specific thallus area (mm\U00B2/dry mg)",
      "mass_to_height" ~ "Mass:height (dry mg/cm)",
      "aspect_ratio" ~ "Thallus length:width",
      "chlA" ~ "chlorophyll A concentration",
      "frond_dmc" ~ "Frond dry matter content",
      "total_dmc" ~ "Thallus dry matter content"
    )
  )) %>% 
  mutate(boxplot = map2(
    data, units,
    ~ ggplot(data = .x,
             aes(x = reorder(x = scientific_name, 
                             X = value,
                             FUN = median,
                             na.rm = TRUE),
                 y = value)) +
      geom_boxplot(aes(fill = sp_code),
                   outliers = FALSE) +
      geom_point(alpha = 0.2,
                 position = position_jitter(width = 0.1,
                                            seed = 666),
                 shape = 21) +
      scale_fill_manual(values = algae_spcode_colors) +
      facet_grid(cols = vars(pigment_type), 
                 scales = "free_x") +
      scale_x_discrete(labels = scales::label_wrap(10)) +
      labs(title = .y,
           x = "Scientific name",
           y = .y) +
      basic_theme
  ))

# FvFm
pluck(boxplots, 4, 1)

# thickness
pluck(boxplots, 4, 2)

# volume
pluck(boxplots, 4, 3)

# SA:P
pluck(boxplots, 4, 4)

# SA:V
pluck(boxplots, 4, 5)

# max height
pluck(boxplots, 4, 6)

# STA
pluck(boxplots, 4, 7)

# mass:height
pluck(boxplots, 4, 8)

# aspect ratio
pluck(boxplots, 4, 9)

# chlorophyll A
pluck(boxplots, 4, 10)

# frond dry matter content
pluck(boxplots, 4, 11)

# thallus dry matter content
pluck(boxplots, 4, 12)

# ⟞ b. saving outputs -----------------------------------------------------

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
  select(specimen_ID, scientific_name,
         maximum_height, mass_to_height, sav_mean, thickness_mm_mean, 
         frond_dmc_mean, sta_mean, sav_mean, sap_mean, fvfm_mean, aspect_ratio_mean, 
         frond_length_mean, frond_width_mean, total_wet, total_dry,
         total_volume, chlA_mean, total_dmc) %>% 
  pivot_longer(maximum_height:total_dmc,
               names_to = "trait",
               values_to = "value") %>% 
  nest(.by = "trait",
       data = everything()) %>% 
  mutate(trait = map(
    trait,
    ~ case_match(
      .x,
      "maximum_height" ~ "max_height",
      "mass_to_height" ~ "mass_to_height",
      "sav_mean" ~ "sav",
      "thickness_mm_mean" ~ "thickness",
      "frond_dmc_mean" ~ "frond_dmc",
      "sta_mean" ~ "sta",
      "sap_mean" ~ "sap",
      "fvfm_mean" ~ "fvfm",
      "aspect_ratio_mean" ~ "aspect_ratio",
      "frond_length_mean" ~ "thallus_length",
      "frond_width_mean" ~ "thallus_width",
      "total_wet" ~ "weight_wet",
      "total_dry" ~ "weight_dry",
      "total_volume" ~ "volume",
      "chlA_mean" ~ "chlA",
      "total_dmc" ~ "total_dmc"
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
  mutate(distribution_plot = pmap(
    list(x = data, y = trait, z = bins),
    function(x, y, z) ggplot(data = x,
                             aes(x = value)) +
      geom_histogram(bins = z,
                     fill = "lightgrey",
                     color = "black") +
      labs(title = paste0(y, " (no transformation)")) +
      basic_theme
  )) %>% 
  mutate(distribution_plot_log = pmap(
    list(x = data, y = trait, z = bins),
    function(x, y, z) ggplot(data = x,
                             aes(x = log(value))) +
      geom_histogram(bins = z,
                     fill = "lightgrey",
                     color = "black") +
      labs(title = paste0(y, " (natural log transform)")) +
      basic_theme
  )) %>% 
  mutate(distribution_plot_sqrt = pmap(
    list(x = data, y = trait, z = bins),
    function(x, y, z) ggplot(data = x,
                             aes(x = sqrt(value))) +
      geom_histogram(bins = z,
                     fill = "lightgrey",
                     color = "black") +
      labs(title = paste0(y, " (square root transform)")) +
      basic_theme
  )) %>% 
  mutate(qq = pmap(
    list(x = data, y = trait),
    function(x, y) ggplot(data = x,
                          aes(sample = value)) +
      geom_qq_line(color = "darkgrey") +
      geom_qq(shape = 21) +
      labs(title = paste0(y, " (no transformation)")) +
      basic_theme
  )) %>% 
  mutate(qq_log = pmap(
    list(x = data, y = trait),
    function(x, y) ggplot(data = x,
                          aes(sample = log(value))) +
      geom_qq_line(color = "darkgrey") +
      geom_qq(shape = 21) +
      labs(title = paste0(y, " (natural log transform)")) +
      basic_theme
  )) %>% 
  mutate(qq_sqrt = pmap(
    list(x = data, y = trait),
    function(x, y) ggplot(data = x,
                          aes(sample = sqrt(value))) +
      geom_qq_line(color = "darkgrey") +
      geom_qq(shape = 21) +
      labs(title = paste0(y, " (square root transform)")) +
      basic_theme
  ))

# maximum height
pluck(distributions, 5, 1) # right skewed
pluck(distributions, 8, 1)
pluck(distributions, 6, 1) # log transform looks normal (ish)
pluck(distributions, 9, 1)
pluck(distributions, 7, 1)
pluck(distributions, 10, 1)

# mass to height
pluck(distributions, 5, 2) # right skewed
pluck(distributions, 8, 2)
pluck(distributions, 6, 2) # log transform looks normal
pluck(distributions, 9, 2)
pluck(distributions, 7, 2)
pluck(distributions, 10, 2)

# surface area to volume
pluck(distributions, 5, 3) # normal ish
pluck(distributions, 8, 3)
pluck(distributions, 6, 3) # log transform looks normal
pluck(distributions, 9, 3)
pluck(distributions, 7, 3)
pluck(distributions, 10, 3)

# thickness
pluck(distributions, 5, 4) # right skewed
pluck(distributions, 8, 4)
pluck(distributions, 6, 4) # log transform looks normal
pluck(distributions, 9, 4)
pluck(distributions, 7, 4)
pluck(distributions, 10, 4)

# frond DMC
pluck(distributions, 5, 5) # potentially bimodal
pluck(distributions, 8, 5)
pluck(distributions, 6, 5) 
pluck(distributions, 9, 5)
pluck(distributions, 7, 5)
pluck(distributions, 10, 5)
# no transformations look good

# STA
pluck(distributions, 5, 6) # right skewed
pluck(distributions, 8, 6)
pluck(distributions, 6, 6) 
pluck(distributions, 9, 6) # log transform is better?
pluck(distributions, 7, 6)
pluck(distributions, 10, 6)

# SA:P
pluck(distributions, 5, 7) # right skewed
pluck(distributions, 8, 7)
pluck(distributions, 6, 7) 
pluck(distributions, 9, 7) # log transform is better?
pluck(distributions, 7, 7)
pluck(distributions, 10, 7)

# FvFm
pluck(distributions, 5, 8) # bimodal
pluck(distributions, 8, 8)
pluck(distributions, 6, 8) 
pluck(distributions, 9, 8)
pluck(distributions, 7, 8)
pluck(distributions, 10, 8)
# no transformations look good

# aspect ratio
pluck(distributions, 5, 9) # right skewed
pluck(distributions, 8, 9)
pluck(distributions, 6, 9) 
pluck(distributions, 9, 9) # log transformation looks better?
pluck(distributions, 7, 9)
pluck(distributions, 10, 9)

# frond length
pluck(distributions, 5, 10) # right skewed
pluck(distributions, 8, 10)
pluck(distributions, 6, 10) 
pluck(distributions, 9, 10) # log transformation looks better
pluck(distributions, 7, 10)
pluck(distributions, 10, 10)

# frond width
pluck(distributions, 5, 11) # right skewed
pluck(distributions, 8, 11)
pluck(distributions, 6, 11) 
pluck(distributions, 9, 11)
pluck(distributions, 7, 11)
pluck(distributions, 10, 11)
# no transformations look good

# total wet mass
pluck(distributions, 5, 12) # right skewed
pluck(distributions, 8, 12)
pluck(distributions, 6, 12) # log transform looks normal
pluck(distributions, 9, 12)
pluck(distributions, 7, 12)
pluck(distributions, 10, 12)

# total dry mass
pluck(distributions, 5, 13) # right skewed
pluck(distributions, 8, 13)
pluck(distributions, 6, 13) # log transform looks normal
pluck(distributions, 9, 13)
pluck(distributions, 7, 13)
pluck(distributions, 10, 13)

# total DMC
pluck(distributions, 5, 16)
pluck(distributions, 8, 16)
pluck(distributions, 6, 16) # log transform looks normal
pluck(distributions, 9, 16)
pluck(distributions, 7, 16)
pluck(distributions, 10, 16)

# traits to transform (potentially)
# max height, mass to height, SAV, thickness, STA, SAP, aspect ratio, 
# frond length, wet mass, dry mass, total DMC

# ⟞ b. saving outputs -----------------------------------------------------

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

