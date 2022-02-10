# Functions to process raw DHS datafiles

# Categorise severe illness free-text/labels into a standard set
severe_illness_categorise <- function(data, symptom_var){
  illness_index <- grepl(symptom_var, names(data))
  data <- data[,illness_index]
  names <- get_labs(data)
  names[grepl("weakness", names, ignore.case = TRUE)] <- "Weakness"
  names[grepl("jaundice", names, ignore.case = TRUE)] <- "Jaundice"
  names[grepl("yellow", names, ignore.case = TRUE)] <- "Jaundice"
  names[grepl("bleeding", names, ignore.case = TRUE)] <- "Bleeding"
  names[grepl("breathing", names, ignore.case = TRUE)] <- "Impaired_breathing"
  names[grepl("breaths", names, ignore.case = TRUE)] <- "Impaired_breathing"
  names[grepl("respiratory", names, ignore.case = TRUE)] <- "Impaired_breathing"
  names[grepl("no", names, ignore.case = TRUE)] <- "None"
  names[grepl("urine", names, ignore.case = TRUE)] <- "Dark_urine"
  names[grepl("consciousness", names, ignore.case = TRUE)] <- "Impaired_consciousness"
  names[grepl("conciousness", names, ignore.case = TRUE)] <- "Impaired_consciousness"
  names[grepl("conscuiousness", names, ignore.case = TRUE)] <- "Impaired_consciousness"
  names[grepl("conscience", names, ignore.case = TRUE)] <- "Impaired_consciousness"
  names[grepl("consiousness", names, ignore.case = TRUE)] <- "Impaired_consciousness"
  names[grepl("convulsion", names, ignore.case = TRUE)] <- "Seizures"
  names[grepl("seizures", names, ignore.case = TRUE)] <- "Seizures"
  names[grepl("epileptic", names, ignore.case = TRUE)] <- "Seizures"
  names[grepl("heart", names, ignore.case = TRUE)] <- "Heart_problems"
  names[grepl("cardiological", names, ignore.case = TRUE)] <- "Heart_problems"
  names[grepl("drink", names, ignore.case = TRUE)] <- "Not_drinking"
  names[grepl("appetite", names, ignore.case = TRUE)] <- "Not_eating"
  names[grepl("inability to eat", names, ignore.case = TRUE)] <- "Not_eating"
  names[grepl("vomit", names, ignore.case = TRUE)] <- "Vomiting"
  names[grepl("prostration", names, ignore.case = TRUE)] <- "Prostration"
  names[grepl("pale", names, ignore.case = TRUE)] <- "Pale_or_cold"
  
  colnames(data) <- names
  
  new <- names[!names %in% c("Weakness", "Jaundice", "Bleeding",
                             "Impaired_breathing", "None", "Dark_urine",
                             "Impaired_consciousness", "Seizures",
                             "Heart_problems", "Not_drinking", "Vomiting",
                             "Prostration", "Pale_or_cold", "Not_eating")]
  if(length(new) > 0){
    message("New category identified: ", new, ". ")
    message("Consider adding to severe_illness_categorise() function")
  }
  
  return(data)
}

# Process household member data
## pr_data: PR DHS file
process_pr <- function(pr_data, symptom_var){
  
  # Severe illness symptoms 
  ## Not mutually exclusive so coded as individual variables)
  ## Slightly different codes across surveys so search by label
  sev_illness <- severe_illness_categorise(pr_data, symptom_var)
  
  pr_data <- pr_data %>%
    rename(
      country_code = hv000,
      cluster = hv001,
      household = hv002,
      child_line_number = hvidx,
      hh_wealth = hv270,
      car_or_truck = hv212,
      motorcycle_or_scooter = hv211,
      mothers_education = hc61,
      sex = hv104,
      year = hv007,
      slept_last_night = hv103,
      hh_selected_hemoglobin = hv042,
      rdt = hml35,
      microscopy = hml32,
      hh_sample_weight = hv005,
      age_in_months = hc1,
      hemoglobin = hc55
    ) %>% 
    mutate(
      # Deal with missing values
      hc56 = ifelse(hc56 == 999, NA, hc56),
      hb = hc56 / 10,
      anemia_level = ifelse(hb < 5, "severe", "non_severe")) %>%
    dplyr::select(country_code, cluster, household, child_line_number,
                  year, hh_selected_hemoglobin, hh_sample_weight,
                  hh_wealth, car_or_truck, motorcycle_or_scooter,  mothers_education, 
                  slept_last_night, sex, age_in_months,
                  rdt, microscopy,
                  hemoglobin,hb, anemia_level
    ) %>%
    bind_cols(sev_illness) %>%
    # replace codes with factor labels
    mutate_if(is.labelled, as_factor) %>%
    # Hyper variables
    mutate(
      # Adjust for RDT results labelled differently (e.g. "falciparum positive")
      rdt = forcats::fct_relabel(rdt, function(x){
        ifelse(grepl("positive", x), "positive", x)
      }),
      age_in_months = as.integer(age_in_months),
      cluster = as.integer(cluster),
      age_group = cut(age_in_months, breaks = c(0, 6, 12, 24, 36, 48, 60), right = FALSE),
      motorised_vehicle = ifelse(car_or_truck == "yes" | motorcycle_or_scooter == "yes", "yes", "no"),
      year = as.factor(year),
      hh_wealth = recode(hh_wealth, lowest = "poorest", second = "poorer", fourth = "richer", highest = "richest"))
  
  pr_data$severe_malaria <- NA
  if("rdt" %in% names(pr_data) & all(c("anemia_level", "Impaired_consciousness", "Seizures", "Impaired_breathing") %in% names(pr_data))){
    message("Severe malaria specified")
    pr_data <- pr_data %>%
      mutate(severe_malaria = ifelse(rdt == "positive" & (anemia_level == "severe" | Impaired_consciousness == "yes" | Seizures == "yes" | Impaired_breathing == "yes"), 1, 0))
  }
  
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
      fever = h22,
      blood_taken = h47,
      no_fever_treatment = h32y,
      fever_treatment_location = h46a,
      convenience_hf_location = s1045,
      distance_to_health_facility = v467d) %>%
    mutate(antimalarial = 
             case_when(ml13e == 1 ~ "ACT",
                       ml13a == 1 ~ "SP/Fansidar",
                       ml13b == 1 ~ "Chloroquine",
                       ml13c == 1 ~ "Amodiaquine",
                       ml13d == 1 ~ "Quinine pills",
                       ml13da == 1 ~ "Quinine injection or IV",
                       ml13aa == 1 ~ "Rectal artesunate",
                       ml13ab == 1 ~ "Artesunate injection or IV",
                       ml13h == 1 ~ "Other",
                       TRUE ~ "None"),
           act = ifelse(antimalarial == "ACT", "yes", "no")) %>%
    dplyr::select(country_code, cluster, household, child_line_number, alive,
                  fever, blood_taken, no_fever_treatment, fever_treatment_location, antimalarial, act,
                  convenience_hf_location, distance_to_health_facility) %>%
    # replace codes with factor labels
    mutate_if(is.labelled, as_factor) %>%
    mutate(child_line_number = as.integer(child_line_number),
           cluster = as.integer(cluster),)
  
  # For a small number of surveys there are a small number (<5) duplicated records. Remove them.
  kr_data <- kr_data[!duplicated(kr_data[,c("cluster", "country_code", "household", "child_line_number")]),]
  
  return(kr_data)
}

# Process Geographic data
## shape: Shape file for country
## prevalence: PfPR prevalence raster
## travel_time_foot: Travel time to nearest HF/hospital on foot raster
## travel_time_vehicle: Travel time to nearest HF/hospital by vehicle raster
## Buffer radius: 5km buffer around geolocations for summarising rasters suggested on pg 30 of "Guidelines on the Use of DHS GPS Data"
process_ge <- function(ge_data, shape, prevalence, travel_time_foot, travel_time_vehicle, hf_location, buffer_radius = 5000){
  
  # Mask to country boundary (speeds up extract process a lot)
  prevalence_sub <- raster::mask(raster::crop(prevalence, extent(shape)), shape)
  travel_time_foot_sub <- raster::mask(raster::crop(travel_time_foot, extent(shape)), shape)
  travel_time_vehicle_sub <- raster::mask(raster::crop(travel_time_vehicle, extent(shape)), shape)
  hf_location_sub <- dplyr::filter(hf_location, ISO == unique(shape$ISO))
  ge_data <- sf::st_transform(ge_data, sf::st_crs(hf_location_sub))
  
  ge_data <- ge_data  %>%
    dplyr::rename(
      cluster = DHSCLUST,
      lat = LATNUM,
      long = LONGNUM,
      ur = URBAN_RURA
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
      prevalence_group = cut(prevalence, seq(0, 1, 0.1), right = FALSE),
      # Travel time in hours
      travel_time_foot = raster::extract(x = travel_time_foot_sub, y = ., buffer = buffer_radius, fun = mean) / 60,
      travel_time_vehicle = raster::extract(x = travel_time_vehicle_sub, y = ., buffer = buffer_radius, fun = mean) / 60
    ) 
  
  nearest <- st_nearest_feature(ge_data, hf_location_sub)
  ge_data$distance <- as.numeric(st_distance(ge_data, hf_location_sub[nearest,], by_element = TRUE) / 1000)
  
  ge_data <- ge_data%>%
    dplyr::select(
      admin1, cluster, lat, long, ur, prevalence, prevalence_group, travel_time_foot, travel_time_vehicle, distance
    ) %>%
    # Remove geometry as not used for statistical models
    st_drop_geometry()
}


dhs_filter_select <- function(data, ...){
  data %>%
    filter(
      #alive == "yes",
      age_in_months %in% 6:59,
      slept_last_night == "yes",
      # fever == "yes",
      rdt %in% "positive",
      #hemoglobin == "measured"
    ) %>%
    dplyr::select(...)
}