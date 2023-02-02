# SBC LTER understory algal functional traits

## last updated 2022-01-31

### File structure

#### code

⊣ `00a_libraries.R`: run this first!  
⊣ `00b_google-sheets.R`: run this if any trait data has been entered on google sheet (last update: 2023-02-02)  

- `00_set-up.R`: set up code, which includes libraries, data, aesthetics  
- `01_cleaning-and-summarizing.R`: basic cleaning, sources `00_set-up.R`  
- `02_prelim-visualizations.Rmd`: boxplots, basically  
- `03_prelim-ordinations.Rmd`: very messy PCoA and PCA  
- `04_benthic-data.Rmd`: looking at timeseries of benthic surveys

#### data

- `00-coarse_traits.csv`: old coarse traits data  
- `algae_all.csv`: all algal species from biomass data set with coarse traits  
- `spp_names.csv`: all species with common and scientific names  
- photosynthetic efficiency data:  
  - `2021-06-30-MOHK-v2-cleaned.csv`  
  - `2021-07-09-MOHK-v2-cleaned.csv`  
  - `2021-07-19-IVEE-v2-cleaned.csv`  
  - `2021-07-20-CARP-v2-cleaned.csv`  
  - `2021-07-21-BULL-v2-cleaned.csv`  
  - `2021-07-22-AQUE-v2-cleaned.csv`  

#### figures

Lots of these, so little time to describe...

### Knitted outputs

- [benthic data](https://an-bui.github.io/algae-traits/code/04_benthic-data.html)  
- [preliminary visualizations](https://an-bui.github.io/algae-traits/code/02_prelim-visualizations.html)

