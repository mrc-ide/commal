# Prepare data for model fitting #

# Load packages
library(dplyr)

### Format DHS data ############################################################
dhs_data_raw <- readRDS("ignore/dhs/processed_data/processed_dhs.RDS")#%>%
  #filter(country %in% c("Kenya", "Tanzania", "Uganda"))
country_levels <- unique(dhs_data_raw$country)

dhs_sma <- dhs_data_raw  %>%
  filter(
    rdt %in% c("negative", "positive"),
    anemia_level %in% c("non_severe", "severe"),
    !is.na(prevalence)) %>%
  mutate(
    # Severe malaria anaemia
    sma = ifelse(anemia_level == "severe" & rdt == "positive", 1, 0)
  ) %>%
  rename(pfpr = prevalence) %>%
  dplyr::select(iso, country, cluster, survey_year, pfpr, rdt, act, sma, distance) %>%
  mutate(countryn = as.numeric(factor(country, levels = country_levels)))


dhs_act_coverage <- dhs_data_raw %>%
  filter(
    rdt == "positive",
    act %in% c("no", "yes")) %>%
  group_by(country) %>%
  summarise(act = mean(act == "yes"))

saveRDS(dhs_sma, "ignore/prob_hosp/dhs_sma.rds")
saveRDS(dhs_act_coverage, "ignore/prob_hosp/dhs_act_coverage.rds")
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
  left_join(dhs_act_coverage, by = "country") %>%
  mutate(countryn = as.numeric(factor(country, levels = country_levels)))

saveRDS(paton_data, "ignore/prob_hosp/paton_inferred.rds")
################################################################################