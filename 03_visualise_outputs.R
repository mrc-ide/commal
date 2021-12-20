library(drjacoby)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

source("R/prevalence_incidence.R")
source("R/odds_probability.R")
source("R/paton_fits.R")
source("R/model_functions.R")

# mcmc <- readRDS("ignore/prob_hosp/mcmc_fits/mcmc.rds")
dhs_sma <- readRDS("ignore/prob_hosp/dhs_sma.rds") %>%
  select(pfpr, sma_microscopy, distance, country)
paton <- readRDS("ignore/prob_hosp/paton_inferred.rds") %>%
  select(pfpr, sma, distance, act, py, country)

### Wrangle parameters #########################################################
samples <- sample_chains(mcmc, 100)

global_parameters <- samples %>%
  select(-contains("ccc_"), -contains("hosp"))

country_parameters <- samples %>%
  select(contains("ccc_"), sample) %>%
  pivot_longer(-sample, names_to = "country", values_to = "country_capacity", names_prefix = "ccc_") 

hospital_parameters <- samples %>%
  select(contains("hosp_"), sample) %>%
  pivot_longer(-sample, names_to = "country", values_to = "hosp", names_prefix = "hosp_") 

parameters <- global_parameters %>%
  left_join(country_parameters, by = "sample") %>%
  left_join(hospital_parameters, by = c("country", "sample")) %>%
  select(country, sample, everything())
################################################################################

### Prepare and summarise DHS data #############################################
dhs_var_summary <- dhs_sma %>%
  bind_rows() %>%
  group_by(country) %>%
  summarise(mean_distance = mean(distance),
            max_distance = max(distance),
            mean_pfpr = mean(pfpr),
            max_pfpr = max(pfpr))

dhs_pfpr_summary <- dhs_sma %>%
  bind_rows() %>%
  group_by(country) %>%
  mutate(pfpr_g = cut_number(pfpr, 5)) %>%
  group_by(country, pfpr_g) %>%
  summarise(
    pfprl = quantile(pfpr, 0.025),
    pfpru = quantile(pfpr, 0.975),
    pfpr = mean(pfpr),
    smal = binom::binom.exact(sum(sma_microscopy), n())$lower,
    smau = binom::binom.exact(sum(sma_microscopy), n())$upper,
    sma_prevalence = mean(sma_microscopy))

dhs_distance_summary <- dhs_sma %>%
  bind_rows() %>%
  group_by(country) %>%
  mutate(distance_g = cut_number(distance, 5)) %>%
  group_by(country, distance_g) %>%
  summarise(
    distancel = quantile(distance, 0.025),
    distanceu = quantile(distance, 0.975),
    distance = mean(distance),
    smal = binom::binom.exact(sum(sma_microscopy), n())$lower,
    smau = binom::binom.exact(sum(sma_microscopy), n())$upper,
    sma_prevalence = mean(sma_microscopy))

paton_var_summary <- paton %>%
  bind_rows() %>%
  group_by(country) %>%
  summarise(mean_distance = weighted.mean(distance, py),
            max_distance = max(distance),
            mean_pfpr = weighted.mean(pfpr, py),
            max_pfpr = max(pfpr),
            mean_act = weighted.mean(act, py),
            max_act = max(act))

paton_summary <- paton %>%
  bind_rows() %>%
  mutate(hospital = 1000 * (sma / py))

################################################################################

### Prediction data.frames #####################################################
pfpr_df <- data.frame(pfpr = seq(0, 0.8, 0.01))
distance_df <- data.frame(distance = seq(0.1, 50, 0.1))

dhs_pfpr_prediction <- parameters %>%
  left_join(pfpr_df, by = character()) %>%
  # add mean distance for each country
  left_join(select(dhs_var_summary, country, mean_distance), by = "country")%>%
  rename(distance = mean_distance) %>%
  # Filter PfPr outside of observed range
  left_join(select(dhs_var_summary, country, max_pfpr), by = "country") %>%
  filter(pfpr <= max_pfpr) %>%
  mutate(sma_prevalence = gompertz(pfpr = pfpr,
                                   distance = distance,
                                   global_capacity = global_capacity, 
                                   country_capacity = country_capacity,
                                   distance_beta = distance_beta,
                                   pfpr_beta = pfpr_beta,
                                   shift = shift))

dhs_pfpr_median_prediction <- dhs_pfpr_prediction %>%
  group_by(country, pfpr) %>%
  summarise(sma_prevalence = median(sma_prevalence))



dhs_distance_prediction <- parameters %>%
  left_join(distance_df, by = character()) %>%
  # add mean pfpr for each country
  left_join(select(dhs_var_summary, country, mean_pfpr), by = "country") %>%
  rename(pfpr = mean_pfpr) %>%
  # Filter distance outside of observed range
  left_join(select(dhs_var_summary, country, max_distance), by = "country") %>%
  filter(distance <= max_distance) %>%
  mutate(sma_prevalence = gompertz(pfpr = pfpr,
                                   distance = distance,
                                   global_capacity = global_capacity, 
                                   country_capacity = country_capacity,
                                   distance_beta = distance_beta,
                                   pfpr_beta = pfpr_beta,
                                   shift = shift))

dhs_distance_median_prediction <- dhs_distance_prediction %>%
  group_by(country, distance) %>%
  summarise(sma_prevalence = median(sma_prevalence))

pfpr_prediction <- parameters %>%
  mutate(country_capacity = 0,
         act = 0.5) %>%
  left_join(pfpr_df, by = character()) %>%
  mutate(distance = mean(dhs_var_summary$mean_distance)) %>%
  mutate(sma_prevalence = gompertz(pfpr = pfpr,
                                   distance = distance,
                                   global_capacity = global_capacity, 
                                   country_capacity = country_capacity,
                                   distance_beta = distance_beta,
                                   pfpr_beta = pfpr_beta,
                                   shift = shift))

pfpr_median_prediction <- pfpr_prediction %>%
  group_by(pfpr) %>%
  summarise(sma_prevalence = median(sma_prevalence))


distance_prediction <- parameters %>%
  mutate(country_capacity = 0) %>%
  left_join(distance_df, by = character()) %>%
  mutate(pfpr = mean(dhs_var_summary$mean_pfpr)) %>%
  mutate(sma_prevalence = gompertz(pfpr = pfpr,
                                   distance = distance,
                                   global_capacity = global_capacity, 
                                   country_capacity = country_capacity,
                                   distance_beta = distance_beta,
                                   pfpr_beta = pfpr_beta,
                                   shift = shift))

distance_median_prediction <- distance_prediction %>%
  group_by(distance) %>%
  summarise(sma_prevalence = median(sma_prevalence))



paton_prediction <- parameters %>%
  filter(country %in% paton_var_summary$country) %>%
  left_join(pfpr_df, by = character()) %>%
  left_join(select(paton_var_summary, country, mean_distance, mean_act), by = "country") %>%
  rename(distance = mean_distance,
         act = mean_act,
    ) %>%
  mutate(sma_prevalence = gompertz(pfpr = pfpr,
                                   distance = distance,
                                   global_capacity = global_capacity, 
                                   country_capacity = country_capacity,
                                   distance_beta = distance_beta,
                                   pfpr_beta = pfpr_beta,
                                   shift = shift),
         adjusted_sma_prevalence = sma_prev_age_standardise(sma_prevalence),
         duration = mean_duration(dur_die = dur_die,
                                  dur_recover = dur_recover,
                                  cfr = cfr),
         community = inc1(prevalence = adjusted_sma_prevalence, recovery_rate = 1 / duration),
         hospital = community / hosp,
         p_hosp = 1 - p(hosp))

paton_median_prediction <- paton_prediction %>%
  group_by(country, pfpr) %>%
  summarise(community = median(community),
            hospital = median(hospital),
            p_hosp = median(p_hosp))

hospital_comunity <- paton_prediction %>%
  select(country, pfpr, sample, community, hospital) %>%
  pivot_longer(c(community, hospital), names_to = "location", values_to = "inc")

hospital_community_median <- hospital_comunity %>%
  group_by(country, pfpr, location) %>%
  summarise(inc = median(inc))

paton_model_prediction <- 
  data.frame(country = c("Kenya", "Tanzania", "Uganda")) %>%
  left_join(pfpr_df, by = character()) %>%
  mutate(hospital = paton_sma(pfpr))
  
################################################################################

### DHS fit plots ##############################################################
dhs_pfpr_fits <- ggplot() +
  geom_line(data = dhs_pfpr_prediction, aes(x = pfpr, y = sma_prevalence, group = sample), alpha = 0.1, col = "#00798c") +
  geom_line(data = dhs_pfpr_median_prediction, aes(x = pfpr, y = sma_prevalence), col = "#edae49", size = 1) +
  geom_linerange(data = dhs_pfpr_summary, aes(x = pfpr, ymin = smal, ymax = smau), col = "#2e4057") +
  geom_linerange(data = dhs_pfpr_summary, aes(y = sma_prevalence, xmin = pfprl, xmax = pfpru), col = "#2e4057") +
  geom_point(data = dhs_pfpr_summary, aes(x = pfpr, y = sma_prevalence), col = "#2e4057") +
  xlab(expression(~italic(Pf)~Pr[2-10])) +
  ylab(expression(SMAPr[0.5-5])) + 
  facet_wrap(~ country, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA))

dhs_distance_fits <- ggplot() +
  geom_line(data = dhs_distance_prediction, aes(x = distance, y = sma_prevalence, group = sample), alpha = 0.1, col = "#00798c") +
  geom_line(data = dhs_distance_median_prediction, aes(x = distance, y = sma_prevalence), col = "#edae49", size = 1) +
  geom_linerange(data = dhs_distance_summary, aes(x = distance, ymin = smal, ymax = smau), col = "#2e4057") +
  geom_linerange(data = dhs_distance_summary, aes(y = sma_prevalence, xmin = distancel, xmax = distanceu), col = "#2e4057") +
  geom_point(data = dhs_distance_summary, aes(x = distance, y = sma_prevalence), col = "#2e4057") +
  xlab("Distance") +
  ylab(expression(SMAPr[0.5-5])) + 
  facet_wrap(~ country, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA))
################################################################################

### Global plots of SMA predictions ############################################
global_pfpr_plot <- ggplot() +
  geom_line(data = pfpr_prediction, aes(x = pfpr, y = sma_prevalence, group = sample), alpha = 0.1, col = "#00798c") +
  geom_line(data = pfpr_median_prediction, aes(x = pfpr, y = sma_prevalence), col = "#edae49", size = 1) +
  xlab(expression(~italic(Pf)~Pr[2-10])) +
  ylab(expression(SMAPr[0.5-5])) + 
  theme_bw()

global_distance_plot <- ggplot() +
  geom_line(data = distance_prediction, aes(x = distance, y = sma_prevalence, group = sample), alpha = 0.1, col = "#00798c") +
  geom_line(data = distance_median_prediction, aes(x = distance, y = sma_prevalence), col = "#edae49", size = 1) +
  xlab("Distance") +
  ylab(expression(SMAPr[0.5-5])) + 
  theme_bw()
################################################################################

### Paton fit plots ############################################################
paton_fit <- ggplot() +
  geom_line(data = paton_prediction, aes(x = pfpr, y = hospital, group = sample), alpha = 0.1, col = "#00798c") +
  geom_line(data = paton_median_prediction, aes(x = pfpr, y = hospital), col = "#edae49", size = 1) +
  geom_line(data = paton_model_prediction, aes(x = pfpr, y = hospital), col = "#d1495b") +
  geom_point(data = paton_summary, aes(x = pfpr, y = hospital), col = "#2e4057") +
  scale_linetype_manual(values = c(1, 2, 1), labels = c("Median", "95% CrI", "draw"), name = "") +
  scale_colour_manual(values = c("red", "darkred", "black"), labels = c("Median", "95% CrI", "draw"), name = "") +
  scale_alpha(guide = "none") + 
  xlab(expression(~italic(Pf)~Pr[2-10])) +
  ylab("Annual hospitalised incidence per 1000 children") +
  facet_wrap(~ country) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA))
################################################################################

### Hospitalisation plot #######################################################
ggplot() +
  geom_line(data = hospital_comunity, aes(x = pfpr, y = inc, col = location, group = interaction(sample, location)), alpha = 0.1) +
  geom_line(data = hospital_community_median, aes(x = pfpr, y = inc, col = location), size = 1) +
  scale_color_manual(values = c("#00798c", "8d96a3")) +
  scale_linetype_manual(values = c(1, 2, 1), labels = c("Median", "95% CrI", "draw"), name = "") +
  scale_alpha(guide = "none") + 
  ylab("Annual incidence per 1000 children") +
  facet_wrap(~country, scales = "free_y") +
  theme_bw() +
  scale_y_log10()

ggplot(paton_prediction, aes(x = p_hosp)) +
  geom_histogram(bins = 30, fill = "darkblue") +
  facet_wrap(~country, scales = "free_y") +
  theme_bw()
################################################################################

### Random effects #############################################################
ggplot(parameters, aes(x = country, fill = country, y = country_capacity)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_boxplot(coef = 100) +
  xlab("") + 
  ylab("") +
  scale_fill_discrete(guide = "none") +
  coord_flip() +
  theme_bw()
################################################################################


### Parameter plots ############################################################
p <- c("global_capacity", "pfpr_beta", "distance_beta", "shift", "group_sd",
       "cfr", "dur_recover", "dur_die", "dur_tx")

# Global parameters
pp <- plot_par(mcmc, p, display = FALSE)
trace_plots <- patchwork::wrap_plots(lapply(pp, function(x) x$trace + theme(legend.position = "none")), ncol = length(p))
hist_plots <- patchwork::wrap_plots(lapply(pp, function(x) x$hist), ncol = length(p))
acf_plots <- patchwork::wrap_plots(lapply(pp, function(x) x$acf), ncol = length(p))
parameter_plots <- trace_plots / hist_plots / acf_plots

# Correlations between global pars
cor <- apply(combn(p, 2), 2, function(x){
  plot_cor(mcmc, x[1], x[2])
})
correlation_plots <- patchwork::wrap_plots(cor) + plot_layout(guides = "collect")
################################################################################