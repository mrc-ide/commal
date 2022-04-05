### PPC - posterior predictive checks against DHS and Paton data ################

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)

source("R/model_functions.R")
source("R/cascade.R")
# Load data and add PPC grouping variables
dhs_sma <- readRDS("ignore/prob_hosp/dhs_sma.rds") %>%
  mutate(country = ifelse(country == "Congo Democratic Republic", "DRC", country)) %>%
  select(country, pfpr, sma) %>%
  mutate(pfprg = cut(pfpr, c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 1), include.lowest = T, labels = c("0-1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5+")))

paton <- readRDS("ignore/prob_hosp/paton_inferred.rds") %>%
  rename(sma = sma_n_diamond) %>%
  select(pfpr, sma, py, distance, country, site, year_start)

# Load fit
parameters <- readRDS("ignore/prob_hosp/mcmc_fits/parameters.rds") %>%
  mutate(country = ifelse(country == "Congo Democratic Republic", "DRC", country)) 

# Prediction - DHS
predict <- parameters %>%
  slice_sample(n = 200) %>%
  select(sample, country, global_capacity, country_capacity, shift, pfpr_beta) %>%
  left_join(dhs_sma, by = "country") %>%
  mutate(predict = pmap_dbl(select(., formalArgs(gompertz)), gompertz)) %>%
  filter(!is.na(predict))

pd1 <- predict %>%
  summarise(ml = mean(sma))

ppc1 <- ggplot() +
  geom_histogram(data = predict, aes(x = predict), fill = "darkblue", bins = 30) +
  geom_vline(data = pd1, aes(xintercept = ml), col = "red") +
  xlab("SMA prevalence") +
  theme_bw()

pd2 <- predict %>%
  group_by(country) %>%
  summarise(ml = mean(sma))

ppc2 <- ggplot() +
  geom_histogram(data = predict, aes(x = predict), fill = "darkblue", bins = 30) +
  geom_vline(data = pd2, aes(xintercept = ml), col = "red") +
  xlab("SMA prevalence") +
  theme_bw() +
  facet_wrap(~country, ncol = 6, scales = "free") +
  theme(strip.background = element_rect(fill = NA))

pd3 <- predict %>%
  group_by(pfprg) %>%
  summarise(ml = mean(sma))

ppc3 <- ggplot() +
  geom_histogram(data = predict, aes(x = predict), fill = "darkblue", bins = 30) +
  geom_vline(data = pd3, aes(xintercept = ml), col = "red") +
  xlab("SMA prevalence") +
  theme_bw() +
  facet_wrap(~pfprg, scales = "free") +
  theme(strip.background = element_rect(fill = NA))

# Prediction - Paton
predict2 <- parameters %>%
  group_by(country) %>%
  #slice_sample(n = 500) %>%
  ungroup() %>%
  filter(country %in% unique(paton$country)) %>%
  left_join(paton, by = "country")  %>%
  mutate(sma_prevalence = pmap_dbl(select(., formalArgs(gompertz)), gompertz),
         predict = cascade(sma_prevalence = sma_prevalence,
                           pfpr = pfpr,
                            chronic = chronic,
                            dur = dur,
                            py = 1000,
                            hosp = hosp,
                            distance_beta = distance_beta,
                            distance = distance)) %>%
  filter(!is.na(predict)) %>%
  mutate(siteyear = paste0(site, year_start))

pd6 <- predict2 %>%
  summarise(ml = mean(1000 * sma / py))

ppc6 <- ggplot() +
  geom_histogram(data = predict2, aes(x = predict), fill = "darkblue", bins = 30) +
  geom_vline(data = pd6, aes(xintercept = ml), col = "red") +
  xlab("Hospitalised incidence per 1000") +
  theme_bw()

pd7 <- predict2 %>%
  group_by(country) %>%
  summarise(ml = mean(1000 * sma / py))

ppc7 <- ggplot() +
  geom_histogram(data = predict2, aes(x = predict), fill = "darkblue", bins = 30) +
  geom_vline(data = pd7, aes(xintercept = ml), col = "red") +
  xlab("Hospitalised incidence per 1000") +
  theme_bw() +
  facet_wrap(~country, scales = "free") +
  theme(strip.background = element_rect(fill = NA))

pd8 <- predict2 %>%
  group_by(siteyear) %>%
  summarise(ml = mean(1000 * sma / py))

ppc8 <- ggplot() +
  geom_histogram(data = predict2, aes(x = predict), fill = "darkblue", bins = 30) +
  geom_vline(data = pd8, aes(xintercept = ml), col = "red") +
  xlab("Hospitalised incidence per 1000") +
  theme_bw() +
  facet_wrap(~siteyear, scales = "free") +
  theme(strip.background = element_rect(fill = NA))

ggsave("figures_tables/ppc/ppc1.png", ppc1)
ggsave("figures_tables/ppc/ppc2.png", ppc2, height = 10, width = 12)
ggsave("figures_tables/ppc/ppc3.png", ppc3)
ggsave("figures_tables/ppc/ppc6.png", ppc6)
ggsave("figures_tables/ppc/ppc7.png", ppc7, height = 3, width = 8)
ggsave("figures_tables/ppc/ppc8.png", ppc8, height = 10, width = 12)
