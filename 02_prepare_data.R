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
    # Non-malarial anemia / Chronic malaria: see page 13 of doi:10.1111/tmi.12313
    chronic_anaemia = ifelse(anemia_level == "severe" & microscopy == "negative", 1, 0),
    # Severe anaemia + malaria + associated symptoms
    symp_sma_microscopy = ifelse(anemia_level == "severe" & microscopy == "positive" & 
                                   (fever == "yes" | Jaundice == "yes" | Pale_or_cold == "yes"), 1, 0)
  ) %>%
  filter(!is.na(symp_sma_microscopy)) %>%
  rename(pfpr = prevalence) %>%
  dplyr::select(iso, country, cluster, survey_year, year, pfpr, microscopy, hb, symp_sma_microscopy, chronic_anaemia) %>%
  mutate(countryn = as.numeric(factor(country, levels = country_levels)))

nrow(dhs_sma)

saveRDS(dhs_sma, "ignore/prob_hosp/dhs_sma.rds")
################################################################################

### Estimate expected number of SMA cases (in lieu of full dataset) ############
paton_data_raw <- read.csv("ignore/paton/paton_table_data.csv")
paton_data_extracted <- read.csv("ignore/paton/paton_sma_adjusted_extract.csv")

# Link to extracted data and convert rates to N
paton_data <- paton_data_raw %>%
  left_join(paton_data_extracted, by = c("site", "year_start", "year_end")) %>%
  drop_na(sma_diamond) %>%
  mutate(
    sma_n_modelled = round((sma_modelled / 1000) * py),
    sma_n_diamond = round((sma_diamond / 1000) * py),
    distance = dist_min + ((dist_max - dist_min) / 2)) %>%
  mutate(countryn = as.numeric(factor(country, levels = country_levels)))

pd <- paton_data  %>%
  mutate(site_date = paste(site, year_start, year_end)) %>%
  select(site_date, pfpr, sma_modelled, sma_diamond) %>%
  mutate(site_date = factor(site_date),
         site_date = forcats::fct_reorder(site_date, sma_modelled)) %>%
  pivot_longer(-c(site_date, pfpr), names_to = "group", values_to = "sma")

ggplot(pd, aes(x = pfpr, y = sma, col = group)) + 
  geom_point() +
  theme_bw()

# Recreate plot S3A from Paton et al:
ggplot(pd, aes(x = site_date, y = sma, col = group)) + 
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

saveRDS(paton_data, "ignore/prob_hosp/paton_inferred.rds")
################################################################################