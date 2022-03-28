### Figure 2 - Hospitalisation fit wrt PfPr ####################################

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)

source("R/model_functions.R")
source("R/prevalence_incidence.R")
source("R/paton_fits.R")
source("R/odds_probability.R")


################################################################################
# Model fit ####################################################################
################################################################################

# Number of uncertainty draws to plot
ndraw <- 200

# Paton data
paton <- readRDS("ignore/prob_hosp/paton_inferred.rds") %>%
  rename(sma = sma_n_modelled) %>%
  select(pfpr, sma, distance, py, country, sma_modelled, sma_diamond)
# Paton model fit
paton_model_prediction <- data.frame(pfpr = seq(0, 0.8, 0.01)) %>%
  mutate(hospital = paton_sma(pfpr)) %>%
  mutate(model = "Paton")

# Load fit
parameters <- readRDS("ignore/prob_hosp/mcmc_fits/parameters.rds") %>%
  filter(country %in% c("Uganda", "Tanzania", "Kenya"))



# Summarise attributed of countries
country_attributes <- paton %>%
  group_by(country) %>%
  summarise(
    distance = mean(distance),
    weight = sum(py)
  )

# Median
median_country_fit <- paton %>%
  left_join(parameters, by = "country") %>%
  select(country, global_capacity, country_capacity, pfpr_beta, shift, chronic, dur, prob_recognise, hosp, distance_beta) %>%
  group_by(country) %>%
  summarise_all(median) %>%
  left_join(country_attributes, by = "country") %>%
  left_join(data.frame(pfpr = seq(0, 0.8, 0.01)), by = character()) %>%
  mutate(mp = pmap_dbl(select(., formalArgs(gompertz)), gompertz),
         malaria_attributable_sma = mp * (1 - chronic),
         symptomatic_sma_prevalence = sma_prev_age_standardise(malaria_attributable_sma),
         community_symptomatic_sma_inc = inc1(prevalence = symptomatic_sma_prevalence, recovery_rate = 1 / dur),
         community_recognised_sma = community_symptomatic_sma_inc * prob_recognise,
         hospital = community_recognised_sma * exp(hosp + distance_beta * distance))

median_global_fit <- median_country_fit %>%
  group_by(pfpr) %>%
  summarise(hospital = weighted.mean(hospital, weight)) %>%
  mutate(model = "Winskill")

# Uncertainty draws
draw_country_fit <- paton %>%
  left_join(parameters, by = "country") %>%
  select(country, global_capacity, country_capacity, pfpr_beta, shift, chronic, dur, prob_recognise, hosp, distance_beta) %>%
  group_by(country) %>%
  slice_sample(n = ndraw) %>%
  mutate(sample = 1:ndraw) %>%
  ungroup() %>%
  left_join(country_attributes, by = "country") %>%
  left_join(data.frame(pfpr = seq(0, 0.8, 0.01)), by = character()) %>%
  mutate(mp = pmap_dbl(select(., formalArgs(gompertz)), gompertz),
         malaria_attributable_sma = mp * (1 - chronic),
         symptomatic_sma_prevalence = sma_prev_age_standardise(malaria_attributable_sma),
         community_symptomatic_sma_inc = inc1(prevalence = symptomatic_sma_prevalence, recovery_rate = 1 / dur),
         community_recognised_sma = community_symptomatic_sma_inc * prob_recognise,
         hospital = community_recognised_sma * exp(hosp + distance_beta * distance))

draw_global_fit <- draw_country_fit %>%
  group_by(pfpr, sample) %>%
  summarise(hospital = weighted.mean(hospital, weight)) %>%
  mutate(model = "Winskill")


# Plotting 
global_plot <- ggplot() +
  geom_point(data = paton, aes(x = pfpr, y = sma_modelled)) +
  geom_line(data = draw_global_fit, aes(x = pfpr, y = hospital, group = sample), alpha = 0.25, col = "#00798c") +
  geom_line(data = median_global_fit, aes(x = pfpr, y = hospital, col = model), col = "#edae49", size = 1) +
  geom_line(data = paton_model_prediction, aes(x = pfpr, y = hospital, col = model), col = "grey20", lty = 2) +
  xlab(expression(~italic(Pf)~Pr[2-10])) +
  ylab("Annual hospitalised incidence per 1000 children") +
  theme_bw() +
  ylim(0, 3.5)

country_plot <- ggplot() +
  geom_point(data = paton, aes(x = pfpr, y = sma_modelled)) +
  geom_line(data = draw_country_fit, aes(x = pfpr, y = hospital, group = sample), alpha = 0.25, col = "#00798c") +
  geom_line(data = median_country_fit, aes(x = pfpr, y = hospital, col = model), col = "#edae49", size = 1) +
  xlab(expression(~italic(Pf)~Pr[2-10])) +
  ylab("Annual hospitalised incidence per 1000 children") +
  theme_bw() +
  facet_wrap(~ country) +
  theme(strip.background = element_rect(fill = NA))
################################################################################
################################################################################
################################################################################

################################################################################
# Prob hospital ################################################################
################################################################################

prob_hosp_pd <- parameters %>%
  mutate(ph = p(exp(hosp)))

# Plotting
prob_hosp_plot <- ggplot(prob_hosp_pd, aes(x = ph)) +
  geom_histogram(binwidth = 0.05, col = "white", fill = "black") + 
  facet_wrap(~ country) +
  xlab("Probability hospital | symptoms recognised") +
  xlim(0, 1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())

# Table
prob_hosp_table <- parameters %>%
  mutate(ph = p(exp(hosp))) %>%
  group_by(country) %>%
  summarise(probability_hospitall = quantile(ph, 0.025),
            probability_hospital = median(ph),
            probability_hospitalu = quantile(ph, 0.975))

write.csv(prob_hosp_table, "figures_tables/probability_hospital.csv", row.names = FALSE)

################################################################################
################################################################################
################################################################################

################################################################################
# Distance #####################################################################
################################################################################

draw_distance_global_fit <- paton %>%
  left_join(parameters, by = "country") %>%
  select(country, hosp, distance_beta) %>%
  group_by(country) %>%
  slice_sample(n = ndraw) %>%
  mutate(sample = 1:ndraw) %>%
  ungroup() %>% 
  left_join(select(country_attributes, country, weight), by = "country") %>%
  left_join(data.frame(distance = seq(0, 100, 0.1)), by = character()) %>%
  mutate(p_hospital = p(exp(hosp + distance_beta * distance))) %>%
  group_by(distance, sample) %>%
  summarise(p_hospital = weighted.mean(p_hospital, weight))

median_distance_global_fit <- paton %>%
  left_join(parameters, by = "country") %>%
  select(country, hosp, distance_beta) %>%
  group_by(country) %>%
  summarise_all(median) %>% 
  left_join(select(country_attributes, country, weight), by = "country") %>%
  left_join(data.frame(distance = seq(0, 100, 0.1)), by = character()) %>%
  mutate(p_hospital = p(exp(hosp + distance_beta * distance))) %>%
  group_by(distance) %>%
  summarise(p_hospital = weighted.mean(p_hospital, weight))

distance_plot <- ggplot() +
  geom_line(data = draw_distance_global_fit, aes(x = distance, y = p_hospital, group = sample), alpha = 0.25, col = "#00798c") +
  geom_line(data = median_distance_global_fit, aes(x = distance, y = p_hospital), col = "#edae49", size = 1) +
  ylab("Probability hospital | symptoms recognised") +
  xlab("Distance") +
  theme_bw()

################################################################################
################################################################################
################################################################################
