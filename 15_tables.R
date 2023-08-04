### Tables #####################################################################

library(dplyr)
library(tidyr)
source("R/odds_probability.R")

options(scipen = 500)
p1 <-  readRDS("ignore/prob_hosp/mcmc_fits/parameters.rds") %>%
  mutate(run = "main")
p2 <-  readRDS("ignore/prob_hosp/mcmc_fits/parameters_sensitivity_fever.rds") %>%
  mutate(run = "fever")
p3 <-  readRDS("ignore/prob_hosp/mcmc_fits/parameters_sensitivity_chronic.rds") %>%
  mutate(run = "chronic")
p4 <-  readRDS("ignore/prob_hosp/mcmc_fits/parameters_sensitivity_dur.rds") %>%
  mutate(run = "dur")
p5 <-  readRDS("ignore/prob_hosp/mcmc_fits/parameters_sensitivity_rdt.rds") %>%
  mutate(run = "rdt")

out <- bind_rows(p1, p2, p3, p4, p5) %>%
  mutate(run = factor(run, levels = c("main", "fever", "rdt", "chronic", "dur"))) 

################################################################################
### Probability hospital #######################################################
################################################################################

prob_hosp <- out %>%
  filter(country %in% c("Uganda", "Tanzania", "Kenya")) %>%
  mutate(ph = p(exp(hosp))) %>%
  group_by(run, country) %>%
  summarise(
    probability_hospitall = quantile(ph, 0.025),
    probability_hospital = median(ph),
    probability_hospitalu = quantile(ph, 0.975))

prob_hosp_table <- prob_hosp %>%
  mutate(formatted = paste0(100 * round(probability_hospital, 2),
                            " (", 100 * round(probability_hospitall, 2),
                            ", ", 100 * round(probability_hospitalu, 2), ")")) %>%
  select(run, country, formatted) %>%
  pivot_wider(country, names_from = run, values_from = formatted)
  
write.csv(prob_hosp_table, "ignore/figures_tables/probability_hospital.csv", row.names = FALSE)
################################################################################
################################################################################
################################################################################

################################################################################
### Parameter summary tables ###################################################
################################################################################
## Global parameters
global_summary <- out %>%
  # We just need one country ad global parameters are replicated for each
  filter(country == "Kenya") %>%
  select(-c(country, chronic, hosp, country_capacity)) %>%
  select(-sample) %>%
  pivot_longer(cols = -run, names_to = "par", values_to = "x") %>%
  group_by(run, par) %>%
  summarise(x = quantile(x, c(0.025, 0.5, 0.975)), q = c(0.05, 0.5, 0.975)) %>%
  ungroup() %>%
  pivot_wider(c(run, par), names_from = q, values_from = x) %>%
  mutate(formatted = paste0(signif(`0.5`, 2),
                            " (", signif(`0.05`, 2),
                            ", ", signif(`0.975`, 2), ")")) %>%
  select(run, par, formatted) %>%
  pivot_wider(par, names_from = run, values_from = formatted)

## 3 Country parameters
country_summary <- out %>%
  filter(country %in% c("Uganda", "Tanzania", "Kenya")) %>%
  select(country, run, chronic, hosp) %>%
  pivot_longer(cols = -c(country, run), names_to = "par", values_to = "x") %>%
  group_by(country, run, par) %>%
  summarise(x = quantile(x, c(0.025, 0.5, 0.975)), q = c(0.05, 0.5, 0.975), na.rm = TRUE) %>%
  ungroup() %>%
  pivot_wider(c(country, run, par), names_from = q, values_from = x) %>%
  mutate(formatted = paste0(signif(`0.5`, 2),
                            " (", signif(`0.05`, 2),
                            ", ", signif(`0.975`, 2), ")")) %>%
  select(country, run, par, formatted) %>%
  pivot_wider(c(country, par), names_from = run, values_from = formatted) %>%
  arrange(desc(par), country)

## All country parameters
re_summary <- out %>%
  select(country, run, country_capacity) %>%
  pivot_longer(cols = -c(country, run), names_to = "par", values_to = "x") %>%
  group_by(country, run, par) %>%
  summarise(x = quantile(x, c(0.025, 0.5, 0.975)), q = c(0.05, 0.5, 0.975), na.rm = TRUE) %>%
  ungroup() %>%
  pivot_wider(c(country, run, par), names_from = q, values_from = x) %>%
  mutate(formatted = paste0(signif(`0.5`, 2),
                            " (", signif(`0.05`, 2),
                            ", ", signif(`0.975`, 2), ")")) %>%
  select(country, run, par, formatted) %>%
  pivot_wider(c(country, par), names_from = run, values_from = formatted)


write.csv(global_summary, "ignore/figures_tables/global_summary.csv", row.names = FALSE)
write.csv(country_summary, "ignore/figures_tables/country_summary.csv", row.names = FALSE)
write.csv(re_summary, "ignore/figures_tables/re_summary.csv", row.names = FALSE)
options(scipen = 0)
################################################################################
################################################################################
################################################################################ 


################################################################################
### Data for source data #######################################################
################################################################################
source_data_prob_hosp <- out %>%
  filter(country %in% c("Uganda", "Tanzania", "Kenya")) %>%
  mutate(ph = p(exp(hosp))) |>
  select(run, country, sample, ph)
write.csv(source_data_prob_hosp, "ignore/figures_tables/source_data/table_1_source_prob_hosp.csv", row.names = FALSE)
################################################################################
################################################################################
################################################################################ 