# This is a script to load in packages

############################################################-
# 1. general use --------------------------------------------
############################################################-

library(tidyverse) # general use
library(here) # file organization

############################################################-
# 2. cleaning -----------------------------------------------
############################################################-

library(janitor) # cleaning names in data frames
library(lubridate) # dealing with dates

############################################################-
# 3. analysis -----------------------------------------------
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
library(nlme) # non-linear mixed effect models (gives significance for terms)
library(MuMIn) # lms
library(glmmTMB) # GLMM
library(emmeans) # effect sizes
library(DHARMa) # plotting residuals from `glmmTMB`
library(performance) # checks of collinearity, amongst other things
library(car) # checking VIR

############################################################-
# 4. visualization ------------------------------------------
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



