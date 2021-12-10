# Wrapper to process country/survey DHS data
## iso: iso3c code
## year: Survey year
## dhs_cache_address: Address of DHS datasets
## shapefile_cache_address: Address of shapefile data
process_country <- function(iso, year, country, symptom_var, dhs_cache_address, shapefile_cache_address){
  
  id_datasets <- read.csv(paste0(dhs_cache_address, "metadata.csv")) %>%
    dplyr::filter(ISO == iso, SurveyYear == year)
  
  pr_data <- readRDS(paste0(dhs_cache_address, "cache/datasets/", id_datasets$household_PR, ".RDS")) %>%
    as_tibble()
  kr_data <- readRDS(paste0(dhs_cache_address, "cache/datasets/", id_datasets$childrens_KR, ".RDS")) %>%
    as_tibble()
  ge_data <- readRDS(paste0(dhs_cache_address, "cache/datasets/", id_datasets$geo_GE, ".RDS")) %>%
    sf::st_as_sf()
  
  # Load country boundary
  shape <- readRDS(shapefile_cache_address) %>%
    filter(ISO == iso)
  # Load rasters
  prevalence <- raster::raster(paste0("ignore/map/getRaster/pfpr_", year, ".tiff"))
  travel_time_foot <- raster::raster(x = "ignore/map/travel_time/2020_walking_only_travel_time_to_healthcare.geotiff")
  travel_time_vehicle <- raster::raster(x = "ignore/map/travel_time/2020_motorized_travel_time_to_healthcare.geotiff")
  hf_location <- readRDS("ignore/map/who-cds-gmp-2019-01-eng.rds")
  
  pr <- process_pr(pr_data, symptom_var)
  kr <- process_kr(kr_data)
  kr <- kr[!duplicated(kr[,c("cluster", "country_code", "household", "child_line_number")]),]
  ge <- process_ge(ge_data, shape, prevalence, travel_time_foot, travel_time_vehicle, hf_location)
  
  data <- ge %>%
    left_join(pr, by = "cluster") %>%
    left_join(kr, by = c("cluster", "country_code", "household", "child_line_number")) %>%
    # Add hypervariable travel time
    mutate(travel_time = ifelse(motorised_vehicle == "yes", travel_time_vehicle, travel_time_foot),
           iso = iso,
           country = country)
  data$survey_year = year
  
  if(nrow(data) > max(nrow(pr), nrow(kr), nrow(ge))){
    stop("Data duplication likely")
  }
  
  return(data)
}