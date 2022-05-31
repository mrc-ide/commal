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
  mutate(sma_prevalence = pmap_dbl(select(., formalArgs(gompertz)), gompertz),
         hospital = cascade(sma_prevalence = sma_prevalence,
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
  mutate(sma_prevalence = pmap_dbl(select(., formalArgs(gompertz)), gompertz),
         hospital = cascade(sma_prevalence = sma_prevalence,
                            pfpr = pfpr,
                            chronic = chronic,
                            dur = dur,
                            py = 1000,
                            hosp = hosp,
                            distance_beta = distance_beta,
                            distance = distance))

median_global_fit <- all_parameters_median %>%
  left_join(data.frame(pfpr = pfpr_range), by = character()) %>%
  mutate(sma_prevalence = pmap_dbl(select(., formalArgs(gompertz)), gompertz),
         hospital = cascade(sma_prevalence = sma_prevalence,
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
  mutate(sma_prevalence = pmap_dbl(select(., formalArgs(gompertz)), gompertz),
         hospital = cascade(sma_prevalence = sma_prevalence,
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
  ylab("Annual hospitalised incidence per 1000 children") +
  theme_bw() +
  ylim(0, 3.5)

ggsave("ignore/figures_tables/Paton_fit.png", global_plot, height = 4, width = 5)

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
################################################################################
################################################################################
################################################################################

################################################################################
# Prob hospital ################################################################
################################################################################

prob_hosp_pd <- country_parameters %>%
  mutate(ph = p(exp(hosp)))

# Plotting
prob_hosp_plot <- function(prob_hosp_pd, title){
  ggplot(prob_hosp_pd, aes(x = ph)) +
  geom_histogram(binwidth = 0.02, col = "white", fill = "black", size = 0.01) + 
  #facet_wrap(~ country) +
  xlab("Hospitalisation probability") +
  ylab(" \n ") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), limits = c(0, 1)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
    ggtitle(title)
}
ph1 <- prob_hosp_plot(filter(prob_hosp_pd, country == "Kenya"), "Kenya")
ph2 <- prob_hosp_plot(filter(prob_hosp_pd, country == "Tanzania"),"Tanzania")
ph3 <- prob_hosp_plot(filter(prob_hosp_pd, country == "Uganda"), "Uganda")

ph4 <- ggplot(prob_hosp_pd, aes(x = country, y = 100 * ph)) +
  geom_boxplot(coef = 100) + 
  xlab("") +
  ylab("% hospitalised") +
  ylim(0, 100) +
  theme_bw()

################################################################################
################################################################################
################################################################################

################################################################################
### Community versus hospital inc ##############################################
################################################################################

community_hosp <- all_parameters_median %>% 
  left_join(data.frame(pfpr = pfpr_range), by = character()) %>%
  mutate(sma_prevalence = pmap_dbl(select(., formalArgs(gompertz)), gompertz),
         ma_sma = malaria_attributable(sma_prevalence, chronic_sa_prevalence = chronic, pfpr = pfpr),
         as_ma_sma = sma_prev_age_standardise(ma_sma),
         Community = prev_to_inc(as_ma_sma, recovery_rate = 1 / dur),
         Hospital = hospitalised(Community, hosp = hosp,
                                 distance_beta = distance_beta,
                                 distance = distance),
         All = Community + Hospital)

spd <- community_hosp %>%
  select(pfpr, Community, Hospital) %>%
  pivot_longer(cols = -pfpr, names_to = "type", values_to = "inc")
         
community_hospital_burden_plot <- ggplot(spd, aes(x = pfpr, y = inc, fill = type)) +
  geom_area() +
  scale_fill_manual(name = "", values = c("#619cff", "#00c19f")) +
  ylab("Annual incidence\nper 1000 children") +
  xlab(expression(~italic(Pf)~Pr[2-10])) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.85),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(0, 0, 0, 0, unit = "cm"),
        legend.title = element_blank()) +
  ggtitle("")

#hospital_summary <- ((ph1 / ph3) | ( ph3 / community_hospital_burden_plot)) + 
 # plot_annotation(tag_levels = "A")

hospital_summary <- ph4 | community_hospital_burden_plot + 
  plot_annotation(tag_levels = "A")

ggsave("ignore/figures_tables/hospital_summary.png", hospital_summary, height = 3, width = 7, scale = 0.8)
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

