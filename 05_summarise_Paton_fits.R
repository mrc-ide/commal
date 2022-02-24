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

paton <- readRDS("ignore/prob_hosp/paton_inferred.rds") %>%
  select(pfpr, sma, distance, py, country)

# Load fit
parameters <- readRDS("ignore/prob_hosp/mcmc_fits/parameters.rds") %>%
  filter(country %in% c("Uganda", "Tanzania", "Kenya"))

compare <- paton %>%
  left_join(parameters, by = "country") %>%
  mutate(sma_prevalence = pmap_dbl(select(., formalArgs(gompertz)), gompertz),
         malaria_attributable_sma = sma_prevalence * (1 - chronic),
         symptomatic_sma_prevalence = sma_prev_age_standardise(sma_prevalence),
         community_symptomatic_sma_inc = inc1(prevalence = symptomatic_sma_prevalence, recovery_rate = 1 / dur, py),
         community_recognised_sma = community_symptomatic_sma_inc * prob_recognise,
         hospital = community_recognised_sma * exp(hosp + distance_beta * distance)) %>%
  group_by(country, pfpr, distance, sma) %>%
  summarise(hospital = median(hospital))

compare_plot <- ggplot(compare, aes(y = hospital, x = sma, col = country)) +
  geom_point() +
  geom_abline() +
  xlab("Observed SMA") +
  ylab("Predicted SMA") +
  theme_bw() +
  xlim(0, 90) + 
  ylim(0, 90) + 
  coord_fixed()

paton_summary <- paton %>%
  group_by(country) %>%
  summarise(py = sum(py),
            distance = median(distance))

paton_draw <- parameters %>%
  group_by(country) %>%
  slice_sample(n = 100) %>%
  ungroup() %>%
  select(-group_sd) %>%
  left_join(paton_summary, by = "country") %>%
  left_join(data.frame(pfpr = seq(0, 0.8, 0.01)), by = character()) %>%
  mutate(sma_prevalence = pmap_dbl(select(., formalArgs(gompertz)), gompertz),
         malaria_attributable_sma = sma_prevalence * (1 - chronic),
         symptomatic_sma_prevalence = sma_prev_age_standardise(sma_prevalence),
         community_symptomatic_sma_inc = inc1(prevalence = symptomatic_sma_prevalence, recovery_rate = 1 / dur),
         community_recognised_sma = community_symptomatic_sma_inc * prob_recognise,
         hospital = community_recognised_sma * exp(hosp + distance_beta * distance))

paton_median <- paton_draw %>%
  group_by(pfpr, country) %>%
  summarise(hospital = median(hospital))

paton_data_summary <- paton %>%
  bind_rows() %>%
  mutate(hospital = 1000 * (sma / py),
         hospitall = 1000 * (qpois(0.025, sma)/ py),
         hospitalu = 1000 * (qpois(0.975, sma)/ py))

paton_model_prediction <- data.frame(pfpr = seq(0, 0.8, 0.01)) %>%
  mutate(hospital = paton_sma(pfpr))

paton_plot <- ggplot() +
  geom_line(data = paton_draw, aes(x = pfpr, y = hospital, group = sample), alpha = 0.2, col = "#00798c") +
  geom_line(data = paton_median, aes(x = pfpr, y = hospital, col = country), size = 1) +
  geom_point(data = paton_data_summary, aes(x = pfpr, y = hospital, col = country)) +
  #geom_linerange(data = paton_data_summary, aes(x = pfpr, ymin = hospitall, ymax = hospitalu, col = country)) +
  xlab(expression(~italic(Pf)~Pr[2-10])) +
  ylab("Annual hospitalised incidence per 1000 children") +
  facet_wrap(~ country) +
  coord_cartesian(ylim = c(0, 5)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA))

paton_combined <- parameters %>%
  select(-group_sd) %>%
  left_join(data.frame(pfpr = seq(0, 0.8, 0.01)), by = character()) %>%
  left_join(paton_summary) %>%
  # Py weighted country-specific variables:
  mutate(sma_prevalence = pmap_dbl(select(., formalArgs(gompertz)), gompertz),
         malaria_attributable_sma = sma_prevalence * (1 - chronic),
         symptomatic_sma_prevalence = sma_prev_age_standardise(sma_prevalence),
         community_symptomatic_sma_inc = inc1(prevalence = symptomatic_sma_prevalence, recovery_rate = 1 / dur),
         community_recognised_sma = community_symptomatic_sma_inc * prob_recognise,
         hospital = community_recognised_sma * exp(hosp + distance_beta * distance))

paton_draws_combined <- paton_combined %>%
  group_by(pfpr, sample) %>%
  summarise(hospital = spatstat.geom::weighted.median(hospital, py))
  
paton_median_combined <- paton_combined %>%
  group_by(pfpr) %>%
  summarise(hospital = spatstat.geom::weighted.median(hospital, py))

paton_plot_combined <- ggplot() +
  geom_line(data = paton_draws_combined, aes(x = pfpr, y = hospital, group = sample), alpha = 0.2, col = "#00798c") +
  geom_line(data = paton_median_combined, aes(x = pfpr, y = hospital), col = "#edae49", size = 1) +
  geom_line(data = paton_model_prediction, aes(x = pfpr, y = hospital), col = "grey20", lty = 2) +
  geom_point(data = paton_data_summary, aes(x = pfpr, y = hospital, col = country)) +
  xlab(expression(~italic(Pf)~Pr[2-10])) +
  ylab("Annual hospitalised incidence per 1000 children") +
  coord_cartesian(ylim = c(0, 3.2)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA))

fig2_supplement <- (paton_plot | compare_plot + theme(legend.position = "none")) +
  plot_layout(widths = c(3, 1), guides = "collect")
ggsave("figures_tables/paton_supplement.png", fig2_supplement, width = 12, height = 4)

ggsave("figures_tables/fig2.png", paton_plot_combined, width = 7, height = 5)


prob_hosp_pd <- parameters %>%
  mutate(ph = p(exp(hosp)))

ggplot(prob_hosp_pd, aes(x = ph, fill = country)) +
  geom_histogram(binwidth = 0.05) + 
  facet_wrap(~ country) +
  theme_bw()

prob_hosp_table <- parameters %>%
  mutate(ph = p(exp(hosp))) %>%
  group_by(country) %>%
  summarise(probability_hospitall = quantile(ph, 0.025),
            probability_hospital = median(ph),
            probability_hospitalu = quantile(ph, 0.975))

write.csv(prob_hosp_table, "figures_tables/probability_hospital.csv", row.names = FALSE)


### Distance trend:

distance_trend_draw <- parameters %>%
  filter(country %in% c("Tanzania", "Uganda", "Kenya")) %>%
  #group_by(country) %>%
  slice_sample(n = 300) %>%
  mutate(sample = 1:n()) %>%
  #ungroup() %>%
  select(sample, country, hosp, distance_beta) %>%
  left_join(data.frame(distance = seq(0, 50, 0.1)), by = character()) %>%
  mutate(p_hospital = p(exp(hosp + distance_beta * distance)))
distance_trend_median <- distance_trend_draw %>%
  group_by(distance) %>%
  summarise(p_hospital = median(p_hospital))  
  
distance_plot <- ggplot() +
  geom_line(data = distance_trend_draw, aes(x = distance, y = p_hospital, group = sample), alpha = 0.2, col = "#00798c") +
  geom_line(data = distance_trend_median, aes(x = distance, y = p_hospital), size = 1) +
  xlab("Distance") +
  ylab("P") +
  #facet_wrap(~ country, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA))