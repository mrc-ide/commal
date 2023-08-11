### Figure 3 - Hospitalisation and incidence ###################################

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

set.seed(123)

################################################################################
# Model draws###################################################################
################################################################################

# Number of uncertainty draws
ndraw <- 100
samples <- sample(1:30000, ndraw)
# PfPr
pfpr_values <- c(0.1, 0.2, 0.3, 0.4, 0.5)

# Summarise attributes of countries
country_attributes <- readRDS("ignore/prob_hosp/paton_inferred.rds") %>%
  group_by(country) %>%
  summarise(
    distance = weighted.mean(distance, py),
    weight = sum(py))

# Pull out draws across our three countries
fits <- c("", "_sensitivity_fever", "_sensitivity_rdt", "_sensitivity_chronic", "_sensitivity_dur")
fit_names <- c("main", "sensitivity_fever", "sensitivity_rdt", "sensitivity_chronic", "sensitivity_dur")
fit_type <- c("Main", "Sensitivity", "Sensitivity", "Sensitivity", "Sensitivity")

country_parameters <- list()
for(i in seq_along(fits)){
  country_parameters[[i]] <- readRDS(paste0("ignore/prob_hosp/mcmc_fits/parameters", fits[i], ".rds")) %>%
    filter(country %in% c("Uganda", "Tanzania", "Kenya"),
           sample %in% samples) %>%
    left_join(country_attributes, by = "country") %>%
    mutate(name = fit_names[i],
           type = fit_type[i]
    )
}
country_parameters <- bind_rows(country_parameters)
################################################################################
################################################################################
################################################################################

################################################################################
### Hospital probability plot ##################################################
################################################################################
prob_hosp_pd <- country_parameters %>%
  mutate(ph = p(exp(hosp)))

# This includes the sensitivity draws
prob_hosp_plot <- ggplot() +
  geom_jitter(data = filter(prob_hosp_pd, type == "Main"),
              aes(x = country, y = 100 * ph, col = country),
              size = 0.05) +
  geom_jitter(data = filter(prob_hosp_pd, type != "Main"),
              aes(x = country, y = 100 * ph, col = country),
              size = 0.05, alpha = 0.1) +
  geom_boxplot(data = filter(prob_hosp_pd, type == "Main"),
               aes(x = country, y = 100 * ph),
               coef = 100, fill = NA, col = "black", size = 0.3) + 
  xlab("") +
  ylab("% hospitalised") +
  ylim(0, 100) +
  scale_colour_manual(values = c("#8d77ff",  "#ff8658", "#0095a2")) +
  theme_bw() +
  guides(colour = "none")
################################################################################
################################################################################
################################################################################

################################################################################
### Summarise hospital and community incidence #################################
################################################################################

# For a summary across all three countries we take the sample-weighted mean
# parameter estimate for each draw of country-specific parameters
all_parameters <- country_parameters %>%
  select(-country) %>%
  summarise(across(everything(), weighted.mean, weights = weights),
            .by = c("sample", "name", "type"))

inc_pd <- all_parameters |>
  cross_join(data.frame(pfpr = pfpr_values)) %>%
  mutate(masa_prevalence = pmap_dbl(select(., formalArgs(gompertz)), gompertz),
         ma_sma = malaria_attributable(masa_prevalence, chronic_sa_prevalence = chronic, pfpr = pfpr),
         as_ma_sma = sma_prev_age_standardise(ma_sma),
         Community = prev_to_inc(as_ma_sma, recovery_rate = 1 / dur),
         Hospital = hospitalised(Community, hosp = hosp,
                                 distance_beta = distance_beta,
                                 distance = distance),
         All = Community + Hospital) |>
  mutate(pfpr = factor(pfpr, levels = c("0.1", "0.2", "0.3", "0.4", "0.5"))) %>%
  select(name, type, pfpr, Community, Hospital) %>%
  pivot_longer(cols = -c("pfpr", "name", "type"), names_to = "location", values_to = "inc")

inc_average_pd <- inc_pd |>
  filter(type == "Main") |>
  summarise(y = median(inc),
            ymin = y,
            ymax = y,
            .by = c("pfpr", "location"))

inc_plot <- ggplot() +
  geom_jitter(data = filter(inc_pd, type == "Main"),
              aes(x = pfpr, y = inc, colour = location, group = location),
              alpha = 0.8, size = 0.05)  +
  geom_jitter(data = filter(inc_pd, type != "Main"),
              aes(x = pfpr, y = inc, colour = location, group = location),
              alpha = 0.1, size = 0.05)  +
  geom_crossbar(data = filter(inc_average_pd, location == "Hospital"),
                aes(x = pfpr, y = y, ymin = ymin, ymax = ymax),
                width = 0.8,
                linewidth = 0.2,
                colour = "#015748")  +
  geom_crossbar(data = filter(inc_average_pd, location == "Community"),
                aes(x = pfpr, y = y, ymin = ymin, ymax = ymax),
                width = 0.8,
                linewidth = 0.2,
                colour = "#ab025d") +
  scale_colour_manual(name = "", values = c("deeppink", "#00c19f")) +
  ylab("Annual incidence\nper 1000 children") +
  xlab(expression(~italic(Pf)~Pr[2-10])) +
  theme_bw() +
  theme(
    legend.background = element_rect(fill = "transparent", colour = "transparent"),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    legend.position = c(0.2, 0.9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.2, "cm"),
    legend.margin = margin(0, 0, 0, 0, unit = "cm"),
    legend.title = element_blank()
  ) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 0.5)))
################################################################################
################################################################################
################################################################################

################################################################################
### Save combined ##############################################################
################################################################################
hospital_summary <- (prob_hosp_plot | inc_plot) + 
  plot_annotation(tag_levels = "A")

ggsave("ignore/figures_tables/hospital_summary.png", hospital_summary,
       height = 3, width = 7, scale = 0.8)
ggsave("ignore/figures_tables/figure_3_hospital_summary.pdf", hospital_summary,
       height = 75, width = 180, units = "mm", scale = 0.9)
################################################################################
################################################################################
################################################################################

################################################################################
### Data for source data #######################################################
################################################################################
source_prob_hosp <- prob_hosp_pd |>
  select(country, ph)

source_inc <- inc_pd |>
  select(pfpr, inc, location)

write.csv(source_prob_hosp, "ignore/figures_tables/source_data/figure_3_source_prob_hosp.csv", row.names = FALSE)
write.csv(source_inc, "ignore/figures_tables/source_data/figure_3_source_model_incidence.csv", row.names = FALSE)
################################################################################
################################################################################
################################################################################
