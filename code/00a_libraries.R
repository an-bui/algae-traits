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

############################################################-
# 3. analysis packages --------------------------------------
############################################################-

# ⊣ a. multivariate analysis --------------------------------

library(vegan) # ordinations etc
library(cluster) # get gower distance on categorical variables
library(FD) # functional diversity
library(pairwiseAdonis) # pairwise comparisons for permanova
library(factoextra) # extract components from multivariate analyses
library(cati) # partitioning species composition/intraspecific variation

# ⊣ b. models -----------------------------------------------

library(lme4) # GLMERs
library(lmerTest) # tests for lm
library(nlme) # non-linear mixed effect models (gives significance for terms)
library(MuMIn) # lms
library(glmmTMB) # GLMM
library(emmeans) # effect sizes
library(DHARMa) # plotting residuals from `glmmTMB`
library(performance) # checks of collinearity, amongst other things
library(car) # checking VIR

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

# ⊣ b. tables -----------------------------------------------

library(gt) # making tables
library(gtsummary) # summary tables for models

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

# ⊣ b. plot aesthetics --------------------------------------

# colors
rhodo_col <- "#781416"
  ochro_col <- "#CC7540"
    chloro_col <- "#6D5A18"
      
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

