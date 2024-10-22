# This is a script to load in packages and create some useful objects

############################################################-
# 1. general use packages -----------------------------------
############################################################-

library(tidyverse) # general use
library(here) # file organization
library(googlesheets4) # pulling data from google sheets

############################################################-
# 2. cleaning packages --------------------------------------
############################################################-

library(janitor) # cleaning names in data frames
library(lubridate) # dealing with dates
library(fuzzyjoin) # fuzzy matching

############################################################-
# 3. analysis packages --------------------------------------
############################################################-

# ⊣ a. multivariate analysis --------------------------------

library(vegan) # ordinations etc
library(cluster) # get gower distance on categorical variables
library(FD) # functional diversity
library(fundiversity) # another functional diversity package
library(pairwiseAdonis) # pairwise comparisons for permanova
library(factoextra) # extract components from multivariate analyses
library(cati) # partitioning species composition/intraspecific variation
library(RVAideMemoire) # pairwise permanova

# ⊣ b. models -----------------------------------------------

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

############################################################-
# 4. visualization packages ---------------------------------
############################################################-

# ⊣ a. plots ------------------------------------------------

library(corrplot) # correlation plots
library(plotly) # for making the species biomass plot
library(patchwork) # putting plots together
library(multcompView) # pairwise comparisons on plots
library(ggeffects) # plot model predictions
library(emmeans) # also plot model predictions
library(ggnewscale) # multiple color scales on ggplot
library(ggridges) # ridge plots
library(GGally) # pairs plots
# source("https://raw.githubusercontent.com/briatte/ggcorr/master/ggcorr.R") # correlation plot
library(ggConvexHull) # convex hulls
library(ggrepel)

# ⊣ b. tables -----------------------------------------------

library(gt) # making tables
library(gtsummary) # summary tables for models
library(flextable) # another making tables option

# ⊣ c. fonts ------------------------------------------------

library(showtext)
font_add_google("Lato", "Lato")
showtext_auto()

############################################################-
# 5. useful objects -----------------------------------------
############################################################-

# ⊣ a. algae vectors ----------------------------------------

algae_common <- c("PH", "PTCA", # Pterygophora californica 
                  "DL", # Desmarestia ligulata
                  "R", # Rhodymenia californica 
                  "CC", # Chondracanthus corymbiferus 
                  "POLA", # Polyneura latissima 
                  "CYOS", # Stephanocystis osmundacea 
                  "FTHR", # Pterosiphonia dendroidea 
                  "CO", # Corallina officinalis var. chilensis 
                  "LX", # Osmundea spectabilis
                  "GS", # Gracilaria spp. 
                  "GR", # Gelidium robustum
                  "BR", # Halymenia spp.
                  "BO", # Bossiella orbigniana 
                  "FB", # Ectocarpaceae spp. 
                  "BF", # Cryptopleura ruprechtiana 
                  "LAFA", # Laminaria farlowii 
                  "CF", # Callophyllis rhynchocarpa 
                  "DP" # Dictyota spp. 
)

# 11 species from Miller et al. 2012 and from conversation with Bob on 2022-01-18
algae_interest <- c("CYOS", # Stephanocystis osmundacea 
                    "LAFA", # Laminaria farlowii 
                    "MAPY", # Macrocystis pyrifera
                    "PH", "PTCA", # Pterygophora californica 
                    "CF", # Callophyllis flabellulata
                    "CC", # Chondracanthus corymbiferus 
                    "GS", # Gracilaria spp. 
                    "POLA", # Polyneura latissima 
                    "FTHR", # Pterosiphonia dendroidea 
                    "R", # Rhodymenia californica 
                    "EGME", # Egregia menziesii
                    "DL" # Desmarestia ligulata
)

# algae list in proposal
# updated 2023-01-27 with new species
algae_proposal <- c("PH", "PTCA", # Pterygophora californica 
                    "BF", # Cryptopleura ruprechtiana                     
                    "CYOS", # Stephanocystis osmundacea                     
                    "DL", # Desmarestia ligulata                   
                    "CC", # Chondracanthus corymbiferus                     
                    "GS", # Gracilaria spp.                     
                    "CO", # Corallina officinalis var. chilensis 
                    "POLA", # Polyneura latissima 
                    "R", # Rhodymenia californica 
                    "GR", # Gelidium robustum
                    "EH", "EGME", # Egregia menziesii
                    "Nandersoniana", # Nienburgia andersoniana
                    "LH", "LAFA", # Laminaria farlowii
                    "DU", # Dictyopteris undulata
                    "DP", # Dictyota
                    "BO" # Bossiella orbigniana
) 

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


algae_proposal_code_factors <- algae_proposal %>% 
  as_factor() %>% 
  fct_relevel(., "BF", "CC", "GS", "CO", "BO", "POLA", "R", "GR", "Nandersoniana",
              "PTCA", "CYOS", "DL", "EGME", "LAFA", "DU", "DP")

algae_proposal_sciname_factors <- as_factor(
  c("Cryptopleura ruprechtiana", 
    "Chondracanthus corymbiferus; Chondracanthus exasperatus", 
    "Gracilaria spp.", 
    "Corallina officinalis", 
    "Bossiella orbigniana", 
    "Polyneura latissima", 
    "Rhodymenia californica", 
    "Gelidium spp.", 
    "Nienburgia andersoniana",
    "Pterygophora californica", 
    "Stephanocystis osmundacea", 
    "Desmarestia ligulata", 
    "Egregia menziesii", 
    "Laminaria farlowii", 
    "Dictyopteris undulata", 
    "Dictyota binghamiae; Dictyota flabellata; Dictyota coriacea"
  )
)

# algae colors

# "#7A2602" "#8A310B" "#9A3D15" "#AA481F" "#BA5429" "#CA6032" "#DA6B3C" "#EA7746" "#FB8350"

algae_colors <- c(
  # reds
  "Cryptopleura ruprechtiana" = "#C68F76", 
  "Chondracanthus corymbiferus; Chondracanthus exasperatus" = "#F5CFBC", 
  "Gracilaria spp." = "#985030", 
  "Corallina officinalis" = "#D6A48D", 
  "Bossiella orbigniana" = "#893B19", 
  "Polyneura latissima" = "#E5B9A4", 
  "Rhodymenia californica" = "#A86547", 
  "Gelidium spp." = "#B77A5F", 
  "Nienburgia andersoniana" = "#7A2602",
  # browns
  "Pterygophora californica" = "#7A6720", 
  "Stephanocystis osmundacea" = "#BA9C30", 
  "Desmarestia ligulata" = "#CFAE35", 
  "Egregia menziesii" = "#8F7825", 
  "Laminaria farlowii" = "#FAD241", 
  "Dictyopteris undulata" = "#A48A2A", 
  "Dictyota binghamiae; Dictyota flabellata; Dictyota coriacea" = "#E4C03B"
)

algae_spcode_colors <- c(
  # reds
  "BF" = "#C68F76", 
  "CC" = "#F5CFBC", 
  "GS" = "#985030", 
  "CO" = "#D6A48D", 
  "BO" = "#893B19", 
  "POLA" = "#E5B9A4", 
  "R" = "#A86547", 
  "GR" = "#B77A5F", 
  "Nandersoniana" = "#7A2602",
  # browns
  "PTCA" = "#7A6720", 
  "CYOS" = "#BA9C30", 
  "DL" = "#CFAE35", 
  "EGME" = "#8F7825", 
  "LAFA" = "#FAD241", 
  "DU" = "#A48A2A", 
  "DP" = "#E4C03B"
)

algae_spcode_full_names <- c(
  # reds
  "BF" = "Cryptopleura\nruprechtiana", 
  "CC" = "Chondracanthus\nspp.", 
  "GS" = "Gracilaria\nspp.", 
  "CO" = "Corallina\nofficinalis", 
  "BO" = "Bossiella\norbigniana", 
  "POLA" = "Polyneura\nlatissima", 
  "R" = "Rhodymenia\ncalifornica", 
  "GR" = "Gelidium\nspp.", 
  "Nandersoniana" = "Nienburgia\nandersoniana",
  # browns
  "PTCA" = "Pterygophora\ncalifornica", 
  "CYOS" = "Stephanocystis\nosmundacea", 
  "DL" = "Desmarestia\nligulata", 
  "EGME" = "Egregia\nmenziesii", 
  "LAFA" = "Laminaria\nfarlowii", 
  "DU" = "Dictyopteris\nundulata", 
  "DP" = "Dictyota\nspp."
)



# colors
rhodo_col <- "#781416"
  ochro_col <- "#CC7540"
    chloro_col <- "#6D5A18"

# rhodophyta: 1 and 2
cluster1 <- "#CC5A17"
cluster2 <- "#8C5332"
# ochrophyta: 3 and 4
cluster3 <- "#16CCCB"
cluster4 <- "#357777"

cat_cluster1 <- "#CC8613"
cat_cluster2 <- "#1274CC"

frond_dmc_col <- "#BE5A47"
total_dmc_col <- "#BE5A47"
thickness_col <- "#BE5A47"
sta_col <- "#BD973D"
sav_col <- "#BD973D"
fvfm_col <- "#b2d8d8"
max_height_col <- "#4D5B75"
sap_col <- "#CC7556"
mass_to_height_col <- "#4D5B75"
aspect_ratio_col <- "#CC7556"
volume_col <- "#BD973D"

trait_color_palette <- c(

  `Frond DMC` = frond_dmc_col,
  `Total DMC` = total_dmc_col,
  `Thickness` = thickness_col,
  `STA` = sta_col,
  `SA:V` = sav_col,
  `Fv/Fm` = fvfm_col,
  

  `Maximum height` = max_height_col,
  `SA:P` = sap_col,
  `Mass:height` = mass_to_height_col,
  `Aspect ratio` = aspect_ratio_col
)
      
    # gradient palette
    gradient_palette <- c("#FFFFFF", "#009BB0")
      
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

