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
source("R/cascade.R")

################################################################################
# Model fit ####################################################################
################################################################################

# Number of uncertainty draws to plot
ndraw <- 100

pfpr_range <- seq(0, 0.7, 0.01)

# Paton data
paton <- readRDS("ignore/prob_hosp/paton_inferred.rds") %>%
  rename(sma = sma_n_modelled) %>%
  select(pfpr, sma, distance, py, country, sma_modelled, sma_diamond)
# Paton model fit
paton_model_prediction <- data.frame(pfpr = pfpr_range) %>%
  mutate(hospital = paton_sma(pfpr)) %>%
  mutate(model = "Paton")
# Summarise attributes of countries
country_attributes <- paton %>%
  group_by(country) %>%
  summarise(
    distance = weighted.mean(distance, py),
    weight = sum(py))

# Load fit
country_parameters <- readRDS("ignore/prob_hosp/mcmc_fits/parameters.rds") %>%
  filter(country %in% c("Uganda", "Tanzania", "Kenya")) %>%
  left_join(country_attributes, by = "country")
country_parameters_median <- country_parameters %>%
  group_by(country) %>% 
  summarise(across(everything(), median))

# For a summary across all three countries we take the sample-weighted mean parameter estimate for
# each draw of country-specific parameters
all_parameters <- country_parameters %>%
  filter(sample <= 10000) %>%
  select(-country) %>%
  group_by(sample) %>%
  summarise(across(everything(), weighted.mean, weights = weights)) %>%
  ungroup()

all_parameters_median <- all_parameters %>%
  summarise(across(everything(), median))

# Median
median_country_fit <- country_parameters_median %>%
  left_join(data.frame(pfpr = pfpr_range), by = character()) %>%
  mutate(masa_prevalence = pmap_dbl(select(., formalArgs(gompertz)), gompertz),
         hospital = cascade(masa_prevalence = masa_prevalence,
                            pfpr = pfpr,
                            chronic = chronic,
                            dur = dur,
                            py = 1000,
                            hosp = hosp,
                            distance_beta = distance_beta,
                            distance = distance))

draw_country_fit <- country_parameters %>%
  group_by(country) %>%
  slice_sample(n = ndraw) %>%
  ungroup() %>%
  left_join(data.frame(pfpr = pfpr_range), by = character()) %>%
  mutate(masa_prevalence = pmap_dbl(select(., formalArgs(gompertz)), gompertz),
         hospital = cascade(masa_prevalence = masa_prevalence,
                            pfpr = pfpr,
                            chronic = chronic,
                            dur = dur,
                            py = 1000,
                            hosp = hosp,
                            distance_beta = distance_beta,
                            distance = distance))

median_global_fit <- all_parameters_median %>%
  left_join(data.frame(pfpr = pfpr_range), by = character()) %>%
  mutate(masa_prevalence = pmap_dbl(select(., formalArgs(gompertz)), gompertz),
         hospital = cascade(masa_prevalence = masa_prevalence,
                            pfpr = pfpr,
                            chronic = chronic,
                            dur = dur,
                            py = 1000,
                            hosp = hosp,
                            distance_beta = distance_beta,
                            distance = distance)) %>%
  mutate(model = "Winskill")

draw_global_fit <- all_parameters %>%
  slice_sample(n = ndraw) %>%
  left_join(data.frame(pfpr = pfpr_range), by = character()) %>%
  mutate(masa_prevalence = pmap_dbl(select(., formalArgs(gompertz)), gompertz),
         hospital = cascade(masa_prevalence = masa_prevalence,
                            pfpr = pfpr,
                            chronic = chronic,
                            dur = dur,
                            py = 1000,
                            hosp = hosp,
                            distance_beta = distance_beta,
                            distance = distance)) %>%
  mutate(model = "Winskill")

# Plotting 
global_plot <- ggplot() +
  geom_point(data = paton, aes(x = pfpr, y = sma_modelled)) +
  geom_line(data = draw_global_fit, aes(x = pfpr, y = hospital, group = sample), alpha = 0.25, col = "#00798c") +
  geom_line(data = median_global_fit, aes(x = pfpr, y = hospital, col = model), col = "#edae49", size = 1) +
  geom_line(data = paton_model_prediction, aes(x = pfpr, y = hospital, col = model), col = "grey20", lty = 2) +
  xlab(expression(~italic(Pf)~Pr[2-10])) +
  ylab("Annual hospitalised incidence \nper 1000 children") +
  theme_bw() +
  ylim(0, 3.6)

ggsave("ignore/figures_tables/Paton_fit.png", global_plot, height = 4, width = 5)
ggsave("ignore/figures_tables/figure_2_paton_fit.pdf", global_plot, height = 70, width = 88, units = "mm", scale = 1.3)


country_plot <- ggplot() +
  geom_point(data = paton, aes(x = pfpr, y = sma_modelled)) +
  geom_line(data = draw_country_fit, aes(x = pfpr, y = hospital, group = sample), alpha = 0.25, col = "#00798c") +
  geom_line(data = median_country_fit, aes(x = pfpr, y = hospital, col = model), col = "#edae49", size = 1) +
  xlab(expression(~italic(Pf)~Pr[2-10])) +
  ylab("Annual hospitalised incidence per 1000 children") +
  theme_bw() +
  facet_wrap(~ country) +
  theme(strip.background = element_rect(fill = NA))

ggsave("ignore/figures_tables/Paton_country_fit.png", country_plot, height = 4, width = 8)

median_global_fit$pfpr[which.min(abs(median_global_fit$hospital - 0.5))]
median_global_fit$pfpr[which.min(abs(median_global_fit$hospital - 1))]
median_global_fit$pfpr[which.min(abs(median_global_fit$hospital - 1.5))]
################################################################################
################################################################################
################################################################################

################################################################################
# Distance #####################################################################
################################################################################

draw_distance_global_fit <- all_parameters %>%
  select(-distance)  %>%
  slice_sample(n = ndraw) %>%
  left_join(data.frame(distance = seq(0, 50, 0.1)), by = character()) %>%
  mutate(p_hospital = p(exp(hosp + distance_beta * distance))) %>%
  group_by(distance, sample) %>%
  summarise(p_hospital = weighted.mean(p_hospital, weight))

median_distance_global_fit <- all_parameters_median %>%
  select(-distance)  %>%
  left_join(data.frame(distance = seq(0, 50, 0.1)), by = character()) %>%
  mutate(p_hospital = p(exp(hosp + distance_beta * distance))) %>%
  group_by(distance, sample) %>%
  summarise(p_hospital = weighted.mean(p_hospital, weight))

distance_plot <- ggplot() +
  geom_line(data = draw_distance_global_fit, aes(x = distance, y = p_hospital, group = sample), alpha = 0.25, col = "#00798c") +
  geom_line(data = median_distance_global_fit, aes(x = distance, y = p_hospital), col = "#edae49", size = 1) +
  ylab("Hospitalisation probability") +
  xlab("Distance") +
  theme_bw()

ggsave("ignore/figures_tables/figS_distance.png", distance_plot, height = 4, width = 4)
################################################################################
################################################################################
################################################################################

################################################################################
### Data for source data #######################################################
################################################################################
source_data_median <- median_global_fit |>
  select(pfpr, hospital) |>
  mutate(
    output = "modelled",
    type = "median"
  )
source_data_draws <- draw_global_fit |>
  select(sample, pfpr, hospital) |>
  mutate(
    output = "modelled",
    type = "draw"
  )

source_data_distance_median <- draw_distance_global_fit |>
  select(-sample) |>
  mutate(
    output = "modelled",
    type = "median"
  )

source_data_distance_draw <- median_distance_global_fit |>
  mutate(
    output = "modelled",
    type = "draw"
  )

write.csv(source_data_median, "ignore/figures_tables/source_data/figure_2_source_model_median.csv", row.names = FALSE)
write.csv(source_data_draws, "ignore/figures_tables/source_data/figure_2_source_model_draws.csv", row.names = FALSE)
write.csv(source_data_distance_median, "ignore/figures_tables/source_data/figure_S1_source_model_median.csv", row.names = FALSE)
write.csv(source_data_distance_draw, "ignore/figures_tables/source_data/figure_S1_source_model_draws.csv", row.names = FALSE)
################################################################################
################################################################################
################################################################################

