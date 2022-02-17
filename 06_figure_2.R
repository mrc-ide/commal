### Figure 2 - Hospitalisation fit wrt PfPr ####################################

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)

source("R/model_functions.R")
source("R/prevalence_incidence.R")
source("R/paton_fits.R")

paton <- readRDS("ignore/prob_hosp/paton_inferred.rds") %>%
  select(pfpr, sma, distance, fever_tx, py, country)

# Load fit
parameters <- readRDS("ignore/prob_hosp/mcmc_fits/parameters.rds") %>%
  filter(country %in% c("Uganda", "Tanzania", "Kenya"))

paton_summary <- paton %>%
  group_by(country) %>%
  summarise(distance = weighted.mean(distance, py),
            fever_tx = weighted.mean(fever_tx, py),
            py = sum(py))

fit_median <- parameters %>%
  select(-sample, -group_sd) %>%
  left_join(paton_summary, by = "country") %>%
  left_join(data.frame(pfpr = seq(0, 0.8, 0.01)), by = character()) %>%
  mutate(sma_prevalence = pmap_dbl(select(., -country, -py, -prob_symptomatic, -dur, -prob_recognise, -hosp), gompertz),
         adjusted_sma_prevalence = sma_prev_age_standardise(sma_prevalence),
         symptomatic_sma_prevalence = adjusted_sma_prevalence * prob_symptomatic,
         community_symptomatic_sma_inc = inc1(prevalence = symptomatic_sma_prevalence, recovery_rate = 1 / dur),
         community_recognised_sma = community_symptomatic_sma_inc * prob_recognise,
         hospital = community_recognised_sma *  rlogit(hosp  -0.19 * distance)) %>%
  group_by(pfpr, country) %>%
  summarise(hospital = median(hospital))

fit_draw <- parameters %>%
  group_by(country) %>%
  slice_sample(n = 100) %>%
  ungroup() %>%
  select(-group_sd) %>%
  left_join(paton_summary, by = "country") %>%
  left_join(data.frame(pfpr = seq(0, 0.8, 0.01)), by = character()) %>%
  mutate(sma_prevalence = pmap_dbl(select(., -prob_symptomatic, -dur, -prob_recognise, -hosp, -sample, -country, -py), gompertz),
         adjusted_sma_prevalence = sma_prev_age_standardise(sma_prevalence),
         symptomatic_sma_prevalence = adjusted_sma_prevalence * prob_symptomatic,
         community_symptomatic_sma_inc = inc1(prevalence = symptomatic_sma_prevalence, recovery_rate = 1 / dur),
         community_recognised_sma = community_symptomatic_sma_inc * prob_recognise,
         hospital = community_recognised_sma * hosp)

paton_data_summary <- paton %>%
  bind_rows() %>%
  mutate(hospital = 1000 * (sma / py))

paton_model_prediction <- data.frame(pfpr = seq(0, 0.8, 0.01)) %>%
  mutate(hospital = paton_sma(pfpr))

ggplot() +
  geom_line(data = fit_draw, aes(x = pfpr, y = hospital, group = sample), alpha = 0.2, col = "#00798c") +
  geom_line(data = fit_median, aes(x = pfpr, y = hospital, col = country), size = 1) +
  geom_line(data = paton_model_prediction, aes(x = pfpr, y = hospital), col = "#d1495b") +
  geom_point(data = paton_data_summary, aes(x = pfpr, y = hospital, col = country)) +
  xlab(expression(~italic(Pf)~Pr[2-10])) +
  ylab("Annual hospitalised incidence per 1000 children") +
  facet_wrap(~ country) +
  coord_cartesian(ylim = c(0, 5)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA))


py <- paton %>%
  group_by(country) %>%
  summarise(py = sum(py))

fit_median_combo <- parameters %>%
  select(-sample, -group_sd) %>%
  left_join(data.frame(pfpr = seq(0, 0.8, 0.01)), by = character()) %>%
  # Py weighted country-specific variables:
  mutate(distance = 10.4, fever_tx = 0.75) %>%#, country_capacity = -0.43, hosp = 0.31) %>%
  mutate(sma_prevalence = pmap_dbl(select(., -country, -prob_symptomatic, -dur, -prob_recognise, -hosp), gompertz),
         adjusted_sma_prevalence = sma_prev_age_standardise(sma_prevalence),
         symptomatic_sma_prevalence = adjusted_sma_prevalence * prob_symptomatic,
         community_symptomatic_sma_inc = inc1(prevalence = symptomatic_sma_prevalence, recovery_rate = 1 / dur),
         community_recognised_sma = community_symptomatic_sma_inc * prob_recognise,
         hospital = community_recognised_sma * hosp) %>%
  group_by(pfpr) %>%
  left_join(py) %>%
  summarise(hospital = spatstat.geom::weighted.median(hospital, py))

fit_combo <- parameters %>%
  select(-group_sd) %>%
  left_join(data.frame(pfpr = seq(0, 0.8, 0.01)), by = character()) %>%
  # Py weighted country-specific variables:
  mutate(distance = 10.4, fever_tx = 0.75) %>%
  mutate(sma_prevalence = pmap_dbl(select(., -sample, -country, -prob_symptomatic, -dur, -prob_recognise, -hosp), gompertz),
         adjusted_sma_prevalence = sma_prev_age_standardise(sma_prevalence),
         symptomatic_sma_prevalence = adjusted_sma_prevalence * prob_symptomatic,
         community_symptomatic_sma_inc = inc1(prevalence = symptomatic_sma_prevalence, recovery_rate = 1 / dur),
         community_recognised_sma = community_symptomatic_sma_inc * prob_recognise,
         hospital = community_recognised_sma * hosp) %>%
  group_by(sample, pfpr) %>%
  left_join(py) %>%
  summarise(hospital = spatstat.geom::weighted.median(hospital, py))

ggplot() +
  geom_line(data = fit_combo, aes(x = pfpr, y = hospital, group = sample), alpha = 0.2, col = "#00798c") +
  geom_line(data = fit_median_combo, aes(x = pfpr, y = hospital), col = "#edae49", size = 1) +
  geom_line(data = paton_model_prediction, aes(x = pfpr, y = hospital), col = "grey40", lty = 2) +
  geom_point(data = paton_data_summary, aes(x = pfpr, y = hospital, col = country)) +
  xlab(expression(~italic(Pf)~Pr[2-10])) +
  ylab("Annual hospitalised incidence per 1000 children") +
  coord_cartesian(ylim = c(0, 3.2)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA))


bf <- mcmc$output %>%
  slice_max(n = 1, order_by = loglikelihood + logprior) %>%
  select(-sample, -group_sd, -chain, -iteration, -phase) %>%
  left_join(data.frame(pfpr = seq(0, 0.8, 0.01)), by = character()) %>%
  # Py weighted country-specific variables:
  mutate(distance = 10.4, fever_tx = 0.75) %>%#, country_capacity = -0.43, hosp = 0.31) %>%
  mutate(sma_prevalence = pmap_dbl(select(., -country, -prob_symptomatic, -dur, -prob_recognise, -hosp), gompertz),
         adjusted_sma_prevalence = sma_prev_age_standardise(sma_prevalence),
         symptomatic_sma_prevalence = adjusted_sma_prevalence * prob_symptomatic,
         community_symptomatic_sma_inc = inc1(prevalence = symptomatic_sma_prevalence, recovery_rate = 1 / dur),
         community_recognised_sma = community_symptomatic_sma_inc * prob_recognise,
         hospital = community_recognised_sma * hosp) %>%
  group_by(pfpr) %>%
  left_join(py) %>%
  summarise(hospital = spatstat.geom::weighted.median(hospital, py))


l1 <- data.frame(sma_prevalence = gompertz(pfpr = seq(0, 0.8, 0.01), distance = 10.4, fever_tx = 0.75,
                                           global_capacity = -4.89, country_capacity = 0.0004, 
                                           distance_beta = 0.154, pfpr_beta = 8.588, tx_beta = -1.31,
                                           shift = 1.535)) %>%
  mutate(adjusted_sma_prevalence = sma_prev_age_standardise(sma_prevalence),
         symptomatic_sma_prevalence = adjusted_sma_prevalence * 0.41,
         community_symptomatic_sma_inc = inc1(prevalence = symptomatic_sma_prevalence, recovery_rate = 1 / 8.78),
         community_recognised_sma = community_symptomatic_sma_inc * 0.276,
         hospital = community_recognised_sma * 0.15,
         ph = hospital / (community_recognised_sma + hospital))

tz <- filter(paton_data_summary, country == "Tanzania")
plot(l1$hospital ~ seq(0, 0.8, 0.01))
points(1000 * tz$sma / tz$py ~ tz$pfpr)
