# Functions to process raw DHS datafiles

# Wrapper to process country/survey DHS data
## iso: iso3c code
## year: Survey year
## dhs_cache_address: Address of DHS datasets
## shapefile_cache_address: Address of shapefile data
process_country <- function(iso, year, country, dhs_cache_address, shapefile_cache_address){
  
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
  prevalence <- raster::raster(paste0("ignore/map/getRaster/pfpr_", min(year, 2019), ".tiff"))
  
  pr <- process_pr(pr_data)
  kr <- process_kr(kr_data)
  ge <- process_ge(ge_data, shape, prevalence)
  
  data <- ge %>%
    left_join(pr, by = "cluster") %>%
    left_join(kr, by = c("cluster", "country_code", "household", "child_line_number")) %>%
    mutate(country = country,
           iso = iso,
           year = year)
  
  if(nrow(data) > max(nrow(pr), nrow(ge))){
    stop("Data duplication likely")
  }
  
  return(data)
}

# Process household member data
## pr_data: PR DHS file
process_pr <- function(pr_data){
  
  pr_data <- pr_data %>%
    rename(
      country_code = hv000,
      cluster = hv001,
      household = hv002,
      child_line_number = hvidx,
      rdt = hml35,
      microscopy = hml32,
      hh_sample_weight = hv005
    ) %>% 
    mutate(
      # Deal with missing values
      hc56 = ifelse(hc56 == 999, NA, hc56),
      hb = hc56 / 10,
      anemia_level = ifelse(hb < 5, "severe", "non_severe")) %>%
    dplyr::select(country_code, cluster, household, child_line_number,
                  hh_sample_weight,
                  rdt, microscopy,
                  hb, anemia_level
    ) %>%
    # replace codes with factor labels
    mutate_if(is.labelled, as_factor) %>%
    # Hyper variables
    mutate(
      # Adjust for RDT results labelled differently (e.g. "falciparum positive")
      rdt = forcats::fct_relabel(rdt, function(x){
        ifelse(grepl("positive", x), "positive", x)
      }),
      cluster = as.integer(cluster))

  
  return(pr_data)
}

# Process Children's data
## kr_data: KR DHS file
process_kr <- function(kr_data){
  possible_missing <- c("s1045", "v467d", "ml13da", "ml13aa", "ml13ab")
  is_missing <- possible_missing[!possible_missing %in% colnames(kr_data)]
  kr_data[,is_missing] <- NA
  
  kr_data <- kr_data %>%
    rename(
      country_code = v000,
      cluster = v001,
      household = v002,
      child_line_number = b16,
      alive = b5,
      fever = h22) %>%
    dplyr::select(country_code, cluster, household, child_line_number, 
                  alive, fever) %>%
    # replace codes with factor labels
    mutate_if(is.labelled, as_factor) %>%
    mutate(child_line_number = as.integer(child_line_number),
           cluster = as.integer(cluster))
  
  # For a small number of surveys there are a small number (<5) duplicated records. Remove them.
  kr_data <- kr_data[!duplicated(kr_data[,c("cluster", "country_code", "household", "child_line_number")]),]
  
  return(kr_data)
}

# Process Geographic data
## shape: Shape file for country
## prevalence: PfPR prevalence raster
## Buffer radius: 5km buffer around geolocations for summarising rasters suggested on pg 30 of "Guidelines on the Use of DHS GPS Data"
process_ge <- function(ge_data, shape, prevalence, buffer_radius = 5000){
  
  # Mask to country boundary (speeds up extract process a lot)
  prevalence_sub <- raster::mask(raster::crop(prevalence, extent(shape)), shape)

  ge_data <- ge_data  %>%
    dplyr::rename(
      cluster = DHSCLUST,
      lat = LATNUM,
      long = LONGNUM
    ) %>%
    # Remove any GIS reference points
    filter(!(round(lat, 6) == 0 & round(long, 6) == 0)) %>%
    # Set country and admin names
    dplyr::mutate(
      admin1 = ifelse(ADM1NAME == "NULL", DHSREGNA, ADM1NAME),
      cluster = as.integer(cluster)
    ) %>%
    # Summarise rasters around clusters
    dplyr::mutate(
      prevalence = raster::extract(x = prevalence_sub, y = ., buffer = buffer_radius, fun = mean),
      prevalence_group = cut(prevalence, seq(0, 1, 0.1), right = FALSE)) 

  ge_data <- ge_data %>%
    dplyr::select(
      admin1, cluster, lat, long, prevalence, prevalence_group) %>%
    # Remove geometry as not used for statistical models
    st_drop_geometry()
  
  return(ge_data)
}
