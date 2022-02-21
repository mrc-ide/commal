### Prepare data for model fitting #############################################

# Load packages
library(dplyr)

### Format DHS data ############################################################
dhs_data_raw <- readRDS("ignore/dhs/processed_data/processed_dhs.RDS")
country_levels <- unique(dhs_data_raw$country)

dhs_sma <- dhs_data_raw  %>%
  filter(
    microscopy %in% c("negative", "positive"),
    anemia_level %in% c("non_severe", "severe"),
    !is.na(prevalence)) %>%
  mutate(
    # Chronic malaria: see page 13 of doi:10.1111/tmi.12313
    chronic_amaemia = ifelse(anemia_level == "severe" & microscopy == "negative", 1, 0),
    # Severe anaemia + malaria + symptoms
    symp_sma_microscopy = ifelse(anemia_level == "severe" & microscopy == "positive" & 
                                   (fever == "yes" | Impaired_consciousness == "yes" | Impaired_breathing == "yes"), 1, 0)
  ) %>%
  filter(!is.na(symp_sma_microscopy)) %>%
  rename(pfpr = prevalence) %>%
  left_join(treatment_access_admin, by = c("country", "admin1")) %>%
  dplyr::select(iso, country, cluster, survey_year, year, pfpr, microscopy, hb, symp_sma_microscopy, chronic_amaemia) %>%
  mutate(countryn = as.numeric(factor(country, levels = country_levels)))


nrow(dhs_sma)

saveRDS(dhs_sma, "ignore/prob_hosp/dhs_sma.rds")
################################################################################

### Estimate expected number of SMA cases (in lieu of full dataset) ############
paton_data_raw <- read.csv("ignore/paton/paton_table_data.csv")
sma_extract <- read.csv("ignore/paton/paton_sma_extract.csv")

nearest_sma <- c()
for(i in 1:nrow(paton_data_raw)){
  nearest_sma[i] <- sma_extract$sma[which.min(abs(paton_data_raw$pfpr[i] - sma_extract$pfpr))]
}

paton_data <- paton_data_raw %>%
  mutate(sma = round((nearest_sma / 1000) * py),
         distance = dist_min + ((dist_max - dist_min) / 2)) %>%
  mutate(countryn = as.numeric(factor(country, levels = country_levels)))

saveRDS(paton_data, "ignore/prob_hosp/paton_inferred.rds")
################################################################################