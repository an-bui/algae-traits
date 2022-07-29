###### 0. source ######

source(here::here("code", "00_set-up.R"))

###### 1. variation in trait values ######

# main question: what is the source of variation in trait values?

###### ⊣ a. fvfm ######

fvfm_nonas <- leaf_traits %>% 
  drop_na(fvfm_mean, sp_code)

# fit a model
fvfm_ind_lme <- lme(log10(fvfm_mean) ~ 1, random = ~1|site/sp_code, data = ind_traits %>% drop_na(fvfm_mean))
fvfm_sub_lme <- lme(log10(fvfm_mean) ~ 1, random = ~1|site/sp_code/specimen_ID/subsample_ID, data = fvfm_nonas)

# variance components
plot(varcomp(fvfm_ind_lme, scale = TRUE))
# site: 69.4%, species: 26.2%, within individuals of same species at same site: 4.4%
plot(varcomp(fvfm_sub_lme, scale = TRUE))
# site: 44.6%, species: 48.8%, specimens: 3.0%, within: 3.7%

###### ⊣ b. STA ######

sta_ind_lme <- lme(log10(sta_mean) ~ 1, random = ~1|site/sp_code, data = ind_traits %>% drop_na(sta_mean))
sta_sub_lme <- lme(log10(sta_mm_mg) ~ 1, random = ~1|site/sp_code/specimen_ID, data = leaf_traits %>% drop_na(sta_mm_mg))

# variance components
plot(varcomp(sta_ind_lme, scale = TRUE))
# species: 75.9%, site: 15.6%, within ind. of same species at same site: 8.5%

plot(varcomp(sta_sub_lme, scale = TRUE))
# species: 69.6, site: 12.8, within a single individual: 12.6%, between individuals: 5.1%

###### ⊣ c. TDMC ######

tdmc_lme <- lme(log10(tdmc_mean) ~ 1, random = ~1|site/sp_code, data = ind_traits %>% drop_na(tdmc_mean))

plot(varcomp(tdmc_lme, scale = TRUE))
# 75.9% by species, 9.2% by site, 6.7% by individual, 8.2% within individual

###### ⊣ d. SA:P ######

sap_lme <- lme(sap_ratio ~ 1, random = ~1|site/sp_code/specimen_ID, data = ind_traits %>% drop_na(sap_ratio))

plot(varcomp(sap_lme, scale = TRUE))
# 62.1% species, 34.9% within individual, 2.2% site, <1% individual

###### ⊣ e. thickness ######

thickness_lme <- lme(thickness_mm_mean ~ 1, random = ~1|site/sp_code/specimen_ID, data = ind_traits %>% drop_na(thickness_mm_mean))

plot(varcomp(thickness_lme, scale = TRUE))
# 58.1% species, 28.4% within individual, 7.7% site, 5.8% individual

###### 2. bootstrapping variance confidence intervals ######

traits <- read.table(here::here("data", "book_data", "chapter6", "Traits.txt"), header=T, stringsAsFactors = T)
commMatrix <- read.table(here::here("data", "book_data", "chapter6","CommMatrix.txt"), header=T, stringsAsFactors = T)

logHeight <- log(traits$Height)
modPart <- lme(logHeight ~ 1, random = ~ 1 | Plot / Species, data = traits, na.action = na.omit)
varcompHeight <- ape::varcomp(modPart, scale = 1)
varcompHeight

nreps <- 1000
varPlot <- varSpecies <- varWithin <- numeric()
for (i in 1:nreps){
  data.aux <- traits[sample(1:nrow(traits), nrow(traits), replace = T), ]
  logHeightAux <- log(data.aux$Height)
  varcomp.aux <- ape::varcomp(lme(logHeightAux ~1 , random = ~1 | Plot / Species, 
                                  data = data.aux, na.action = na.omit), scale = 1)
  varPlot[i] <- varcomp.aux["Plot"]
  varSpecies[i] <- varcomp.aux["Species"]
  varWithin[i] <- varcomp.aux["Within"]
}

quantile(varSpecies, probs=c(0.025, 0.975))

# 1000 iterations
nreps <- 1000

# making holding vectors
varPlot <- varSpecies <- varWithin <- numeric() 

# turn leaf traits into a matrix
ind_traits_mat <- ind_traits %>% 
  column_to_rownames("specimen_ID") %>% 
  select(site, sp_code, fvfm_mean, sta_mean)

# for each iteration
for (i in 1:nreps){
  # resample
  data.aux <- ind_traits_mat[sample(1:nrow(ind_traits_mat), nrow(ind_traits_mat), replace = T), ] 
  # take the log
  logstaAux <- log(data.aux$sta_mean)
  # model with random effects structure
  varcomp.aux <- ape::varcomp(lme(log10(sta_mean) ~ 1, random = ~1|site/sp_code, data = ind_traits %>% drop_na(sta_mean)))
  # site component
  varPlot[i] <- varcomp.aux["site"]
  # species component
  varSpecies[i] <- varcomp.aux["sp_code"]
  # within individual component
  varWithin[i] <- varcomp.aux["Within"]
}
quantile(varSpecies, probs=c(0.025, 0.975))





