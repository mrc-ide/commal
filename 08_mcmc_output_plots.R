
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(drjacoby)

# Load fit
parameters <- readRDS("ignore/prob_hosp/mcmc_fits/parameters.rds") 
mcmc <- readRDS("ignore/prob_hosp/mcmc_fits/mcmc.rds")

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
p <- c(names(select(parameters, -sample, -hosp, -country, -country_capacity)), "hosp_Kenya", "hosp_Tanzania", "hosp_Uganda")

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

### Validation against rtss ####################################################
rtss <- read.csv("ignore/rtss/rtss_trial_data_summary.csv")

#cp <- country_parameters %>%
#  select(country, sample, country_capacity)

rtss_prediction <- parameters %>%
  slice_sample(n = 100) %>%
  select(-country, -hosp, -overdispersion, -group_sd) %>%
  mutate(country_capacity = 0) %>%
  left_join(rtss, by = character()) %>%
  filter(!site == "total") %>%
  mutate(sma_prevalence = gompertz(pfpr = pfpr,
                                   global_capacity = global_capacity, 
                                   country_capacity = country_capacity,
                                   pfpr_beta = pfpr_beta,
                                   shift = shift),
         malaria_attributable_sma = sma_prevalence * (1 - chronic),
         symptomatic_sma_prevalence = sma_prev_age_standardise(malaria_attributable_sma, age_out_lower  = 5, age_out_upper = 17),
         community_symptomatic_sma_inc = inc1(prevalence = symptomatic_sma_prevalence, recovery_rate = 1 / dur, py = py),
         community_recognised_sma = community_symptomatic_sma_inc * prob_recognise) %>%
  group_by(sample) %>%
  summarise(community_symptomatic_sma_inc = sum(community_symptomatic_sma_inc),
            community_recognised_sma = sum(community_recognised_sma)) %>%
  ungroup() %>%
  pivot_longer(-sample, names_to = "measure", values_to = "inc")

ggplot(rtss_prediction, aes(x = measure, y = inc)) +
  geom_boxplot() +
  # Add in RTSS incidence
  geom_hline(yintercept = 54, lty = 2) +
  theme_bw() +
  ylim(0, 500)
################################################################################

