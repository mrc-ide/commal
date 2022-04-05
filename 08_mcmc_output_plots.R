
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
re_plot <- ggplot(parameters, aes(x = country, fill = country, y = country_capacity)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_boxplot(coef = 100) +
  xlab("") + 
  ylab("") +
  scale_fill_discrete(guide = "none") +
  coord_flip() +
  theme_bw()

# Correlation with GPD
re <- parameters %>% 
  select(country, country_capacity, global_capacity) %>%
  group_by(country) %>%
  summarise_all(median) %>%
  mutate(iso2c = countrycode::countrycode(country, "country.name", "iso2c"))

gdp <- WDI::WDI(indicator = "NY.GDP.PCAP.PP.KD", country = re$iso2c, start = 2000, end = 2020) %>%
  select(-country) %>%
  group_by(iso2c) %>%
  filter(!is.na(NY.GDP.PCAP.PP.KD)) %>%
  slice_max(year)

gdp_re <- re %>%
  left_join(gdp, by = "iso2c") %>%
  mutate(upper = rlogit(global_capacity + country_capacity)) %>%
  rename(gdp = NY.GDP.PCAP.PP.KD)

summary(lm(upper ~ log(gdp), data = gdp_re))

re_gdp_plot <- ggplot(gdp_re, aes(x = gdp, y = upper)) +
  geom_smooth(method = 'lm', formula = y~x, se = FALSE, col = "darkred") +
  geom_point() +
  ylim(0, 0.005) +
  xlab("GDP per capita PPP \n (constant 2005 international $)") +
  ylab("Upper SMA prevalence estimate") +
  theme_bw() +
  scale_x_log10()
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

rtss_prediction <- parameters %>%
  slice_sample(n = 100) %>%
  select(-country, -hosp, -overdispersion, -group_sd) %>%
  mutate(country_capacity = 0) %>%
  left_join(rtss, by = character()) %>%
  filter(!site == "total") %>%
  mutate(symptomatic_sma_prevalence = pmap_dbl(select(., formalArgs(gompertz)), gompertz),
         malaria_attributable_sma = malaria_attributable(symptomatic_sma_prevalence, chronic = chronic),
         as_malaria_attributable_sma = sma_prev_age_standardise(malaria_attributable_sma, age_out_lower  = 5, age_out_upper = 17),
         inc = prev_to_inc(as_malaria_attributable_sma, recovery_rate = 1 / dur, py = py),
         recog_inc = recognised(inc, prob_recognise = prob_recognise)) %>%
  group_by(sample) %>%
  summarise(community_symptomatic_sma_inc = sum(inc),
            community_recognised_sma = sum(recog_inc)) %>%
  ungroup() %>%
  pivot_longer(-sample, names_to = "measure", values_to = "inc")

ggplot(rtss_prediction, aes(x = measure, y = inc)) +
  geom_boxplot() +
  # Add in RTSS incidence
  geom_hline(yintercept = 54, lty = 2) +
  theme_bw() +
  ylim(0, 500)
################################################################################

