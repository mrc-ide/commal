### Summarising prob hosp for sensitivity runs #################################

country_parameters_fever <- readRDS("ignore/prob_hosp/mcmc_fits/parameters_sensitivity_fever.rds") %>%
  filter(country %in% c("Uganda", "Tanzania", "Kenya")) %>%
  mutate(ph = p(exp(hosp))) %>%
  group_by(country) %>%
  summarise(
    probability_hospitall_fever = quantile(ph, 0.025),
    probability_hospital_fever = median(ph),
    probability_hospitalu_fever = quantile(ph, 0.975))

country_parameters_chronic <- readRDS("ignore/prob_hosp/mcmc_fits/parameters_sensitivity_chronic.rds") %>%
  filter(country %in% c("Uganda", "Tanzania", "Kenya")) %>%
  mutate(ph = p(exp(hosp))) %>%
  group_by(country) %>%
  summarise(
    probability_hospitall_chronic = quantile(ph, 0.025),
    probability_hospital_chronic = median(ph),
    probability_hospitalu_chronic = quantile(ph, 0.975))

sensitivity_table <- left_join(country_parameters_fever, country_parameters_chronic, by = "country")

write.csv(sensitivity_table, "ignore/figures_tables/probability_hospital_sensitivity.csv", row.names = FALSE)
