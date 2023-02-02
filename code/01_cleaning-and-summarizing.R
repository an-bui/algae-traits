###### 0. source ######

source(here::here("code", "00_set-up.R"))

###### 1. variation in trait values ######

# main question: what is the source of variation in trait values?
# following Messier et al. 2010 and Messier blog post with code: https://juliemessier.org/code/

###### ⊣ a. fvfm ######

fvfm_df <- fvfm_raw %>% 
  rename('fvfm_meas' = '1:Fv/Fm') %>% 
  # fill in '1:Fv/Fm' with 1:Y(II) when there is none
  # only an issue with 20220405-MOHK samples
  mutate(fvfm_meas = case_when(
    !is.na((fvfm_meas)) ~ `1:Y (II)`,
    TRUE ~ fvfm_meas
  )) %>% 
  # take out observations with "-" in the measurement
  filter(fvfm_meas != "-") %>% 
  # make sure that fvfm_meas is numeric
  mutate(fvfm_meas = as.numeric(fvfm_meas)) %>% 
  # only use FvFm, drop the NAs, change the column name
  select(specimen_ID, subsample_ID, fvfm_meas, Time, Date) %>% 
  drop_na(specimen_ID) %>% 
  select(-specimen_ID) %>% 
  left_join(., metadata_sub, by = "subsample_ID") %>% 
  select(site, sp_code, specimen_ID, subsample_ID, fvfm_meas, Time, Date) %>% 
  drop_na() %>% 
  mutate(Date = mdy(Date)) %>% 
  mutate(year = year(Date)) %>% 
  filter(year == 2022)

fvfm_nonas <- leaf_traits %>% 
  drop_na(fvfm_mean, sp_code) %>% 
  filter(year == 2022)

# fit a model
fvfm_ind_lme <- lme(log10(fvfm_mean) ~ 1, random = ~1|site/sp_code, data = ind_traits %>% drop_na(fvfm_mean))
fvfm_sub_lme <- lme(log10(fvfm_mean) ~ 1, random = ~1|site/sp_code/specimen_ID, data = fvfm_nonas)

# variance components
plot(varcomp(fvfm_ind_lme, scale = TRUE))
# species: 89.8%, within individuals: 6.9%, year: 2.9%, site: <1%
plot(varcomp(fvfm_sub_lme, scale = TRUE))
# species: 46.3%, site: 29.3%, within subsample: 18.8%, across individuals: 3.7%, within individuals: 1.93%

###### ⊣ b. STA ######

sta_ind_lme <- lme(log10(sta_mean) ~ 1, random = ~1|site/sp_code, data = ind_traits %>% drop_na(sta_mean) %>% filter(year == 2022))
sta_sub_lme <- lme(log10(sta_mm_mg) ~ 1, random = ~1|site/sp_code/specimen_ID, data = leaf_traits %>% drop_na(sta_mm_mg, sp_code) %>% filter(year == 2022))

# variance components
plot(varcomp(sta_ind_lme, scale = TRUE))
# species: 75.9%, site: 15.6%, within ind. of same species at same site: 8.5%

plot(varcomp(sta_sub_lme, scale = TRUE))
# species: 69.6, site: 12.8, within a single individual: 12.6%, between individuals: 5.1%

###### ⊣ c. TDMC ######

tdmc_ind_lme <- lme(log10(tdmc_mean) ~ 1, random = ~1|site/sp_code, data = ind_traits %>% drop_na(tdmc_mean) %>% filter(year == 2022))

tdmc_sub_lme <- lme(log10(dmc) ~ 1, random = ~1|site/sp_code/specimen_ID, data = leaf_traits %>% drop_na(dmc, sp_code) %>% filter(year == 2022))

plot(varcomp(tdmc_ind_lme, scale = TRUE))
# 78 by species, 8.6 by site, 13.5 within

plot(varcomp(tdmc_sub_lme, scale = TRUE))
# 75.9 by species, 9.2 by site, 6.7 by individual, 8.2 within individual

###### ⊣ d. SA:P ######

sap_ind_lme <- lme(log10(sap_mean) ~ 1, random = ~1|year/site/sp_code/specimen_ID, data = ind_traits %>% drop_na(sap_mean))

sap_sub_lme <- lme(log10(sap_ratio) ~ 1, random = ~1|site/sp_code/specimen_ID, data = leaf_traits %>% drop_na(sap_ratio))

plot(varcomp(sap_ind_lme, scale = TRUE))
# 93.1% species, 2.3% site, 2.9% specimen ID

plot(varcomp(sap_sub_lme, scale = TRUE))
# 82.9% species, 10% site, 2.5% individuals, 4.5% subsample

###### ⊣ e. thickness ######

thickness_df <- thickness %>% 
  pivot_longer(cols = thickness_01:thickness_10, names_to = "measurement_n", values_to = "thickness_mm") %>% 
  left_join(., metadata_sub, by = "subsample_ID") %>% 
  drop_na(thickness_mm) %>% 
  filter(type == "thallus") %>% 
  select(site, sp_code, specimen_ID.x, subsample_ID, thickness_mm) %>% 
  drop_na()

thickness_lme <- lme(thickness_mm_mean ~ 1, random = ~1|site/sp_code/specimen_ID, data = ind_traits %>% drop_na(thickness_mm_mean))

# thickness_sub_lme <- lme(log10(thickness_mm_mean) ~ 1, random = ~1|site/sp_code/specimen_ID/subsample_ID, data = leaf_traits %>% drop_na(thickness_mm_mean))

thickness_sub_lme <- lme(thickness_mm ~ 1, random = ~1|site/sp_code/specimen_ID.x/subsample_ID, data = thickness_df)

plot(varcomp(thickness_lme, scale = TRUE))
# 58.1% species, 28.4% within individual, 7.7% site, 5.8% individual

plot(varcomp(thickness_sub_lme, scale = TRUE))
# 77.5% species, 4.6% individual, 17.9% subsample

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
  filter(year == 2022) %>% 
  select(site, sp_code, fvfm_mean, sta_mean) 

# for each iteration
for (i in 1:nreps){
  # resample
  data.aux <- ind_traits_mat[sample(1:nrow(ind_traits_mat), nrow(ind_traits_mat), replace = T), ] 
  # take the log
  logfvfmAux <- log(data.aux$fvfm_mean)
  # model with random effects structure
  varcomp.aux <- ape::varcomp(lme(logfvfmAux ~1 , random = ~1 | site/sp_code, 
                                  data = data.aux, na.action = na.omit), scale = 1)
  # site component
  varPlot[i] <- varcomp.aux["site"]
  # species component
  varSpecies[i] <- varcomp.aux["sp_code"]
  # within individual component
  varWithin[i] <- varcomp.aux["Within"]
}

varcomp(fvfm_ind_lme, scale = TRUE)

quantile(varSpecies, probs=c(0.025, 0.975))





