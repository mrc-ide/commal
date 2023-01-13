### Prepare data for model fitting #############################################

# Load packages
library(dplyr)

### Format DHS data ############################################################
dhs_data_raw <- readRDS("ignore/dhs/processed_data/processed_dhs.RDS")
 
1 - (mean(is.na(dhs_data_raw$hb) | (is.na(dhs_data_raw$rdt) & is.na(dhs_data_raw$microscopy))))

dhs_masa <- dhs_data_raw  %>%
  filter(
    microscopy %in% c("negative", "positive"),
    !is.na(hb),
    !is.na(prevalence)) %>%
  mutate(
    # Severe anaemia
    sa = ifelse(anemia_level == "severe", 1, 0),
    # Malaria and severe anaemia
    masa  = ifelse(anemia_level == "severe" & microscopy == "positive", 1, 0),
    # Malaria and severe anaemia and fevere
    masa_fever  = ifelse(anemia_level == "severe" & microscopy == "positive" & fever == "yes", 1, 0),
    # Malaria and severe anaemia, RDT defined
    masa_rdt = case_when(anemia_level == "severe" & rdt == "positive" ~ 1,
                        rdt == "negative" ~ 0,
                        anemia_level == "non_severe" ~ 0),
    # Sample weight
    weight = hh_sample_weight / 10^6
  ) %>%
  rename(pfpr = prevalence) %>%
  dplyr::select(iso, country, cluster, year, weight,
                pfpr, microscopy, hb, masa, masa_fever, sa, masa_rdt, rdt)

nrow(dhs_masa)
length(unique(dhs_masa$country))
tapply(dhs_masa$masa, dhs_masa$country, sum)
tapply(dhs_masa$masa, dhs_masa$country, mean)

round(100 * weighted.mean(dhs_masa$sa, dhs_masa$weight), 4)
round(100 * weighted.mean(dhs_masa$microscopy == "positive", dhs_masa$weight), 3)
round(100 * weighted.mean(dhs_masa$masa, dhs_masa$weight), 4)

saveRDS(dhs_masa, "ignore/prob_hosp/dhs_masa.rds")
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