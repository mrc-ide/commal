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
  summarise(distance = mean(distance),
            fever_tx = mean(fever_tx),
            py =sum(py))

fit_median <- parameters %>%
  select(-sample, -group_sd) %>%
  group_by(country) %>%
  summarise_all(.funs = median) %>%
  left_join(paton_summary, by = "country") %>%
  left_join(data.frame(pfpr = seq(0, 0.8, 0.01)), by = character()) %>%
  mutate(sma_prevalence = pmap_dbl(select(., -country, -py, -prob_symptomatic, -dur, -prob_recognise, -hosp), gompertz),
         adjusted_sma_prevalence = sma_prev_age_standardise(sma_prevalence),
         symptomatic_sma_prevalence = adjusted_sma_prevalence * prob_symptomatic,
         community_symptomatic_sma_inc = inc1(prevalence = symptomatic_sma_prevalence, recovery_rate = 1 / dur),
         community_recognised_sma = community_symptomatic_sma_inc * prob_recognise,
         hospital = community_recognised_sma * hosp) %>%
  group_by(pfpr) %>%
  summarise(hospital = median(hospital))

fit_draw <- parameters %>%
  select(-country, -group_sd) %>%
  group_by(sample) %>%
  summarise_all(.funs = median) %>%
  left_join(paton_summary, by = character()) %>%
  left_join(data.frame(pfpr = seq(0, 0.8, 0.01)), by = character()) %>%
  mutate(sma_prevalence = pmap_dbl(select(., -prob_symptomatic, -dur, -prob_recognise, -hosp, -sample), gompertz),
         adjusted_sma_prevalence = sma_prev_age_standardise(sma_prevalence),
         symptomatic_sma_prevalence = adjusted_sma_prevalence * prob_symptomatic,
         community_symptomatic_sma_inc = inc1(prevalence = symptomatic_sma_prevalence, recovery_rate = 1 / dur),
         community_recognised_sma = community_symptomatic_sma_inc * prob_recognise,
         hospital = community_recognised_sma * hosp)

paton_summary <- paton %>%
  bind_rows() %>%
  mutate(hospital = 1000 * (sma / py))

paton_model_prediction <- 
  data.frame(pfpr = seq(0, 0.8, 0.01)) %>%
  mutate(hospital = paton_sma(pfpr))

ggplot() +
  geom_line(data = fit_draw, aes(x = pfpr, y = hospital, group = sample), alpha = 0.2, col = "#00798c") +
  geom_line(data = fit_median, aes(x = pfpr, y = hospital), col = "#edae49", size = 1) +
  geom_line(data = paton_model_prediction, aes(x = pfpr, y = hospital), col = "#d1495b") +
  geom_point(data = paton_summary, aes(x = pfpr, y = hospital, col = country)) +
  xlab(expression(~italic(Pf)~Pr[2-10])) +
  ylab("Annual hospitalised incidence per 1000 children") +
  theme_bw()