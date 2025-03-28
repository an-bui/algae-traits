# This is a script to load in packages and create some useful objects

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------------------------------ 1. packages ------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ⟞ a. general use --------------------------------------------------------

library(tidyverse) # general use
library(here) # file organization
library(googlesheets4) # pulling data from google sheets


# ⟞ b. cleaning -----------------------------------------------------------

library(janitor) # cleaning names in data frames
library(lubridate) # dealing with dates
library(fuzzyjoin) # fuzzy matching


# ⟞ c. multivariate analysis ----------------------------------------------

library(vegan) # ordinations etc
library(cluster) # get gower distance on categorical variables
library(FD) # functional diversity
library(fundiversity) # another functional diversity package
library(pairwiseAdonis) # pairwise comparisons for permanova
library(factoextra) # extract components from multivariate analyses
library(cati) # partitioning species composition/intraspecific variation
library(RVAideMemoire) # pairwise permanova


# ⟞ d. models -------------------------------------------------------------

library(lmerTest) # tests for lm, loads in lme4
library(nlme) # non-linear mixed effect models (gives significance for terms)
library(MuMIn) # lms
library(glmmTMB) # GLMM
library(emmeans) # effect sizes
library(DHARMa) # plotting residuals from `glmmTMB`
library(performance) # checks of collinearity, amongst other things
library(car) # checking VIR
library(specr) # variance decomposition for lmer objects
library(lmodel2) # SMA
library(smatr) # SMA
library(ggpmisc) # visualization of SMA
library(ggExtra) # visualization of SMA
library(broom) # getting tidy model outputs
library(rstatix) # for Dunn test


# ⟞ e. plots --------------------------------------------------------------

library(corrplot) # correlation plots
library(plotly) # for making the species biomass plot
library(patchwork) # putting plots together
library(multcompView) # pairwise comparisons on plots
library(rcompanion) # pairwise comparisons for KW
library(ggeffects) # plot model predictions
library(emmeans) # also plot model predictions
library(ggnewscale) # multiple color scales on ggplot
library(ggridges) # ridge plots
library(GGally) # pairs plots
# source("https://raw.githubusercontent.com/briatte/ggcorr/master/ggcorr.R") # correlation plot
library(ggConvexHull) # convex hulls
library(ggrepel)


# ⟞ f. tables -------------------------------------------------------------

library(gt) # making tables
library(gtsummary) # summary tables for models
library(flextable) # another making tables option


# ⟞ g. fonts --------------------------------------------------------------

library(showtext)
font_add_google("Lato", "Lato")
showtext_auto()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# --------------------------- 2. useful objects ---------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ⟞ a. vectors and colors -------------------------------------------------

# ⟞ ⟞ i. algae ------------------------------------------------------------

algae_colors <- c(
  # BO, CO: articulated calcareous
  "Corallina officinalis" = "#E67932", 
  "Bossiella orbigniana" = "#E49C39", 
  # corticated macrophytes: BF, R
  "Rhodymenia californica" = "#A6B354",
  "Cryptopleura ruprechtiana" = "#77B77D", 
  # corticated foliose: DP
  "Dictyota binghamiae; Dictyota flabellata; Dictyota coriacea" = "#4D8AC6",
  # leathery macrophytes: CC, CYOS, PTCA, LAFA
  "Chondracanthus corymbiferus; Chondracanthus exasperatus" = "#C3ABD1", 
  "Stephanocystis osmundacea" = "#A778B4", 
  "Pterygophora californica" = "#8C4E99", 
  "Laminaria farlowii" = "#6F4C9B"
)

algae_splabel_colors <- c(
  # BO, CO: articulated calcareous
  "Corallina officinalis" = "#E67932", 
  "Bossiella orbigniana" = "#E49C39", 
  # corticated macrophytes: BF, R
  "Rhodymenia californica" = "#A6B354",
  "Cryptopleura ruprechtiana" = "#77B77D", 
  # corticated foliose: DP
  "Dictyota spp." = "#4D8AC6",
  # leathery macrophytes: CC, CYOS, PTCA, LAFA
  "Chondracanthus spp." = "#C3ABD1", 
  "Stephanocystis osmundacea" = "#A778B4", 
  "Pterygophora californica" = "#8C4E99", 
  "Laminaria farlowii" = "#6F4C9B"
)

algae_spcode_colors <- c(
  # BO, CO: articulated calcareous
  "CO" = "#E67932", 
  "BO" = "#E49C39", 
  # corticated macrophytes: BF, R
  "R" = "#A6B354",
  "BF" = "#77B77D", 
  # corticated foliose: DP
  "DP" = "#4D8AC6",
  # leathery macrophytes: CC, CYOS, PTCA, LAFA
  "CC" = "#C3ABD1", 
  "CYOS" = "#A778B4", 
  "PTCA" = "#8C4E99", 
  "LAFA" = "#6F4C9B"
) 

algae_factors <- c(
  # BO, CO
  "Corallina officinalis", 
  "Bossiella orbigniana", 
  # everything else (BF, CC, R, CYOS, DP)
  "Cryptopleura ruprechtiana", 
  "Chondracanthus corymbiferus; Chondracanthus exasperatus", 
  "Rhodymenia californica", 
  "Stephanocystis osmundacea", 
  "Dictyota binghamiae; Dictyota flabellata; Dictyota coriacea",
  # PTCA, LAFA
  "Pterygophora californica", 
  "Laminaria farlowii"
)

algae_splabel_factors <- c(
  "Corallina officinalis", 
  "Bossiella orbigniana", 
  "Rhodymenia californica", 
  "Cryptopleura ruprechtiana", 
  "Dictyota spp.",
  "Chondracanthus spp.", 
  "Stephanocystis osmundacea", 
  "Pterygophora californica", 
  "Laminaria farlowii"
)

algae_spcode_factors <- c(
  "CO", "BO", "R", "BF", "DP", "CC", "CYOS", "PTCA", "LAFA"
)

algae_spcode_full_names <- c(
  # reds
  "BF" = "Cryptopleura ruprechtiana", 
  "CC" = "Chondracanthus spp.", 
  "CO" = "Corallina officinalis", 
  "BO" = "Bossiella orbigniana", 
  "R" = "Rhodymenia californica", 
  # browns
  "PTCA" = "Pterygophora californica", 
  "CYOS" = "Stephanocystis osmundacea", 
  "LAFA" = "Laminaria farlowii", 
  "DP" = "Dictyota spp."
)

# BO and CO are not different from each other = articulated corallines
# DP and BF are not different = foliose ish things
# PTCA and LAFA are not different = stipate browns
# CYOS: bushy stipate brown
# CC: short leathery red
# R: bushy short red
# DP: bushy short brown
groups <- c(
  "BO" = "articulated coralline",
  "CO" = "articulated coralline",
  
  "DP" = "foliose, bushy turf",
  "BF" = "foliose, bushy turf",
  
  "PTCA" = "tall, stipate",
  "LAFA" = "tall, stipate",
  
  "CYOS" = "bushy, stipate",
  
  "CC" = "short, leathery red",
  
  "R" = "bushy, short red"
)


# algae_common <- c("PH", "PTCA", # Pterygophora californica 
#                   "DL", # Desmarestia ligulata
#                   "R", # Rhodymenia californica 
#                   "CC", # Chondracanthus corymbiferus 
#                   "POLA", # Polyneura latissima 
#                   "CYOS", # Stephanocystis osmundacea 
#                   "FTHR", # Pterosiphonia dendroidea 
#                   "CO", # Corallina officinalis var. chilensis 
#                   "LX", # Osmundea spectabilis
#                   "GS", # Gracilaria spp. 
#                   "GR", # Gelidium robustum
#                   "BR", # Halymenia spp.
#                   "BO", # Bossiella orbigniana 
#                   "FB", # Ectocarpaceae spp. 
#                   "BF", # Cryptopleura ruprechtiana 
#                   "LAFA", # Laminaria farlowii 
#                   "CF", # Callophyllis rhynchocarpa 
#                   "DP" # Dictyota spp. 
# )

# 11 species from Miller et al. 2012 and from conversation with Bob on 2022-01-18
# algae_interest <- c("CYOS", # Stephanocystis osmundacea 
#                     "LAFA", # Laminaria farlowii 
#                     "MAPY", # Macrocystis pyrifera
#                     "PH", "PTCA", # Pterygophora californica 
#                     "CF", # Callophyllis flabellulata
#                     "CC", # Chondracanthus corymbiferus 
#                     "GS", # Gracilaria spp. 
#                     "POLA", # Polyneura latissima 
#                     "FTHR", # Pterosiphonia dendroidea 
#                     "R", # Rhodymenia californica 
#                     "EGME", # Egregia menziesii
#                     "DL" # Desmarestia ligulata
# )

# algae list in proposal
# updated 2023-01-27 with new species
# algae_proposal <- c("PH", "PTCA", # Pterygophora californica 
#                     "BF", # Cryptopleura ruprechtiana                     
#                     "CYOS", # Stephanocystis osmundacea                     
#                     "DL", # Desmarestia ligulata                   
#                     "CC", # Chondracanthus corymbiferus                     
#                     "GS", # Gracilaria spp.                     
#                     "CO", # Corallina officinalis var. chilensis 
#                     "POLA", # Polyneura latissima 
#                     "R", # Rhodymenia californica 
#                     "GR", # Gelidium robustum
#                     "EH", "EGME", # Egregia menziesii
#                     "Nandersoniana", # Nienburgia andersoniana
#                     "LH", "LAFA", # Laminaria farlowii
#                     "DU", # Dictyopteris undulata
#                     "DP", # Dictyota
#                     "BO" # Bossiella orbigniana
# ) 

# reds
# "BF", # Cryptopleura ruprechtiana  
# "CC", # Chondracanthus corymbiferus 
# "GS", # Gracilaria spp. 
# "CO", # Corallina officinalis var. chilensis 
# "BO" # Bossiella orbigniana
# "POLA", # Polyneura latissima 
# "R", # Rhodymenia californica 
# "GR", # Gelidium robustum
# "Nandersoniana", # Nienburgia andersoniana

# browns
# "PH", "PTCA", # Pterygophora californica 
# "CYOS", # Stephanocystis osmundacea                     
# "DL", # Desmarestia ligulata                   
# "EH", "EGME", # Egregia menziesii
# "LH", "LAFA", # Laminaria farlowii
# "DU", # Dictyopteris undulata
# "DP", # Dictyota


# algae_proposal_code_factors <- algae_proposal %>% 
#   as_factor() %>% 
#   fct_relevel(., "BF", "CC", "GS", "CO", "BO", "POLA", "R", "GR", "Nandersoniana",
#               "PTCA", "CYOS", "DL", "EGME", "LAFA", "DU", "DP")

# algae_proposal_sciname_factors <- as_factor(
#   c("Cryptopleura ruprechtiana", 
#     "Chondracanthus corymbiferus; Chondracanthus exasperatus", 
#     "Gracilaria spp.", 
#     "Corallina officinalis", 
#     "Bossiella orbigniana", 
#     "Polyneura latissima", 
#     "Rhodymenia californica", 
#     "Gelidium spp.", 
#     "Nienburgia andersoniana",
#     "Pterygophora californica", 
#     "Stephanocystis osmundacea", 
#     "Desmarestia ligulata", 
#     "Egregia menziesii", 
#     "Laminaria farlowii", 
#     "Dictyopteris undulata", 
#     "Dictyota binghamiae; Dictyota flabellata; Dictyota coriacea"
#   )
# )
# blue, green, red, purple


# algae_spcode_colors <- c(
#   # reds
#   "BF" = "#C68F76", 
#   "CC" = "#F5CFBC", 
#   "GS" = "#985030", 
#   "CO" = "#D6A48D", 
#   "BO" = "#893B19", 
#   "POLA" = "#E5B9A4", 
#   "R" = "#A86547", 
#   "GR" = "#B77A5F", 
#   "Nandersoniana" = "#7A2602",
#   # browns
#   "PTCA" = "#7A6720", 
#   "CYOS" = "#BA9C30", 
#   "DL" = "#CFAE35", 
#   "EGME" = "#8F7825", 
#   "LAFA" = "#FAD241", 
#   "DU" = "#A48A2A", 
#   "DP" = "#E4C03B"
# )




# colors
# rhodo_col <- "#781416"
#   ochro_col <- "#CC7540"
#     chloro_col <- "#6D5A18"

# rhodophyta: 1 and 2
# cluster1 <- "#CC5A17"
# cluster2 <- "#8C5332"
# ochrophyta: 3 and 4
# cluster3 <- "#16CCCB"
# cluster4 <- "#357777"
# 
# cat_cluster1 <- "#CC8613"
# cat_cluster2 <- "#1274CC"
# 

# ⟞ ⟞ ii. traits ----------------------------------------------------------

sta_col <- "#BD973D"
sav_col <- "#CC7556"
sap_col <- "#CC7556"
aspect_ratio_col <- "#CC7556"
frond_dmc_col <- "#BE5A47"
total_dmc_col <- "#BE5A47"
thickness_col <- "#BE5A47"
mass_to_height_col <- "#6B6D9F"
max_height_col <- "#4D5B75"
fvfm_col <- "#b2d8d8"
volume_col <- "#CC7556"

height_col <- "#BE5A47"
surface_area_col <- "#BE5A47"
thickness_col <- "#4D5B75"
tradeoff_col <- "#865B5E"

trait_color_palette <- c(
  `Height` = height_col,
  `Thickness` = thickness_col,
  `Surface area` = surface_area_col,
  `Height:wet weight` = tradeoff_col,
  `Dry:wet weight` = tradeoff_col,
  `Height:volume` = tradeoff_col,
  `Surface area:volume` = tradeoff_col,
  `Surface area:dry weight` = tradeoff_col,
  `Surface area:perimeter` = tradeoff_col
)

trait_factor <- c("Height",
                  "Surface area",
                  "Thickness",
                  "Surface area:perimeter",
                  "Surface area:volume",
                  "Surface area:dry weight",
                  "Height:wet weight",
                  "Height:volume",
                  "Dry:wet weight")

trait_abbreviations <- c(
  "H", "SA", "T", "SA:P", "SA:V", "SA:DW", "H:WW", "H:V", "DW:WW"
)

trait_abbreviations <- c(
  "H" = "Height", 
  "SA" = "Surface area", 
  "T" = "Thickness", 
  "SA:P" = "Surface area:perimeter", 
  "SA:V" = "Surface area:volume", 
  "SA:DW" = "Surface area:dry weight", 
  "H:WW" = "Height:wet weight", 
  "H:V" = "Height:volume", 
  "DW:WW" = "Dry:wet weight"
)

trait_colnames_factor <- c(
  "maximum_height",
  "frond_area_scaled",
  "thickness_mm_mean",
  "sap_mean",
  "sav_scaled",
  "sta_scaled",
  "height_ww",
  "height_vol",
  "total_dmc"
)
      
    # gradient palette
    # gradient_palette <- c("#FFFFFF", "#009BB0")
    # 

# ⟞ ⟞ iii. sites ----------------------------------------------------------
      
    # site name vector
    sites <- c("bull", "aque", "ahnd", "napl", "ivee", "golb", 
               "abur", "mohk", "carp", "scdi", "sctw")
    
    sites_proposal <- c("bull", "aque", "napl", "ivee", "mohk", "carp")
    sites_proposal_new <- c("aque", "napl", "ivee", "mohk", "carp")
    
    sites_full <- setNames(c("Bullito (BULL)", 
                             "Arroyo Quemado (AQUE)",
                             "Arroyo Hondo (AHND)",
                             "Naples (NAPL)",
                             "Isla Vista (IVEE)",
                             "Goleta Beach (GOLB)",
                             "Arroyo Burro (ABUR)",
                             "Mohawk (MOHK)",
                             "Carpinteria (CARP)",
                             "Diablo Canyon (SCDI)",
                             "Twin Harbors (SCTW)"), sites)
    
    # full names for sites
    bull_full <- "Bullito (BULL)"
    aque_full <- "Arroyo Quemado (AQUE)"
    ahnd_full <- "Arroyo Hondo (AHND)"
    napl_full <- "Naples (NAPL)"
    ivee_full <- "Isla Vista (IVEE)"
    golb_full <- "Goleta Beach (GOLB)"
    abur_full <- "Arroyo Burro (ABUR)"
    mohk_full <- "Mohawk (MOHK)"
    carp_full <- "Carpinteria (CARP)"
    scdi_full <- "Diablo Canyon (SCDI)"
    sctw_full <- "Twin Harbors (SCTW)"
    


# ⟞ b. functions ----------------------------------------------------------
    
    # function to calculate standard error
    se <- function(x,...){
      sd(x, na.rm = TRUE)/sqrt(length(na.omit(x)))
    }
    
    # function to arrange variance decomp output from specr::icc_specs()
    var_decomp <- function(x) {
      icc_specs(x) %>%
        mutate_if(is.numeric, round, 1) %>% 
        arrange(-percent)
    }

