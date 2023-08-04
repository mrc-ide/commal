### Validation #################################################################

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)

source("R/model_functions.R")
source("R/prevalence_incidence.R")
source("R/paton_fits.R")
source("R/odds_probability.R")
source("R/cascade.R")

# Load fit
parameters <- readRDS("ignore/prob_hosp/mcmc_fits/parameters.rds") 
mcmc <- readRDS("ignore/prob_hosp/mcmc_fits/mcmc.rds")

### Validation against rtss ####################################################
rtss <- read.csv("ignore/rtss/rtss_trial_data_summary.csv")

rtss_prediction <- rtss %>%
  filter(!site == "total") %>%
  select(country, py, pfpr) %>%
  left_join(data.frame(sample = 1:1000), by = character()) %>%
  left_join(unique(select(parameters, sample, global_capacity, shift, pfpr_beta, dur, distance_beta)), by = "sample") %>%
  left_join(select(parameters, country, sample, country_capacity, chronic, hosp), by = c("country", "sample")) %>%
  mutate(chronic = replace_na(chronic, mean(chronic, na.rm = TRUE)), 
         hosp = replace_na(hosp, mean(hosp, na.rm = TRUE)),
         country_capacity = replace_na(country_capacity, 0)) %>%
  mutate(masa_prevalence = pmap_dbl(select(., formalArgs(gompertz)), gompertz),
         malaria_attributable_sma = malaria_attributable(masa_prevalence, chronic = chronic, pfpr = pfpr),
         as_malaria_attributable_sma = sma_prev_age_standardise(malaria_attributable_sma, age_out_lower  = 5, age_out_upper = 17),
         community = prev_to_inc(as_malaria_attributable_sma, recovery_rate = 1 / dur, py = py),
         hospital = hospitalised(community, hosp, distance_beta, 0),
         total = community + hospital) %>%
  group_by(sample) %>%
  summarise(Total = sum(total),
            Hospital = sum(hospital))
 
rtss_pd <- rtss_prediction %>%
  pivot_longer(cols = -sample, names_to = "est", values_to = "y")

rtss_plot <- ggplot(rtss_pd, aes(x = est, y = y)) +
  geom_boxplot(coef = 100) +
  # Add in RTSS incidence
  geom_hline(yintercept = 54, lty = 2, col = "darkred") +
  geom_hline(yintercept = 44, lty = 2, col = "darkblue") +
  ylab("SMA cases") +
  scale_y_log10() +
  xlab("") +
  theme_bw()

ggsave("ignore/figures_tables/rtss_validation.png", rtss_plot, height = 4, width = 4)
################################################################################

################################################################################
### Data for source data #######################################################
################################################################################
write.csv(rtss_pd, "ignore/figures_tables/source_data/figure_S2_source_model.csv", row.names = FALSE)
################################################################################
################################################################################
################################################################################