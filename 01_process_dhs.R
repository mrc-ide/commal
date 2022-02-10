#### Load and process raw DHS survey data ######################################

library(sf)
library(raster)
library(dplyr)
library(tidyr)
library(haven)

# Processing functions
source("R/process_country.R")
source("R/process_dhs.R")
source("R/dhs_helpers.R")

# DHS survey data cache
dropbox <- "C:/Users/pwinskil/Dropbox (SPH Imperial College)/"
dhs_cache <- "DHS_severe_malaria_cache/"

# Choose surveys - we want those with geo-coded clusters and severe symptom variables
surveys <- read.csv(paste0(dropbox, dhs_cache, "metadata.csv")) %>%
  filter(!is.na(symptom_var),
         !is.na(geo_GE),
         SurveyYear < 2020) %>%
  dplyr::select(ISO, CountryName, SurveyYear, symptom_var)

# Processing
data_list <- list()
for(i in 1:nrow(surveys)){
  message("ISO: ", surveys$ISO[i])
  message("Year: ", surveys$SurveyYear[i])
  data_list[[i]] <- process_country(
    iso = surveys$ISO[i], year = surveys$SurveyYear[i], country = surveys$CountryName[i],
    symptom_var = surveys$symptom_var[i],
    dhs_cache_address = paste0(dropbox, dhs_cache),
    shapefile_cache_address = paste0(dropbox, "geoboundaries_cache/sf/admin0.RDS")
  )
}
data <- bind_rows(data_list)
nrow(data)
saveRDS(data, "ignore/dhs/processed_data/processed_dhs.RDS")
