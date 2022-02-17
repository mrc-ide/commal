### PPC - posterior predictive checks against DHS and Paton data ################

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)

source("R/model_functions.R")

# Load data and add PPC grouping variables
dhs_sma <- readRDS("ignore/prob_hosp/dhs_sma.rds") %>%
  mutate(country = ifelse(country == "Congo Democratic Republic", "DRC", country)) %>%
  select(country, pfpr, distance, fever_tx, sma_microscopy) %>%
  mutate(pfprg = cut(pfpr, c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 1), include.lowest = T, labels = c("0-1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5+")),
         distanceg = cut(distance, c(0, 2, 4, 6, 8, 10, 1000), include.lowest = T, labels = c("0-2", "2-4", "4-6", "6-8", "8-10", "10+")),
         txg = cut(fever_tx, c(0, 0.4, 0.6, 0.8, 1), include_lowest = T, labels = c("<0.4", "0.4-0.6", "0.60-0.8", "0.8+")))

paton <- readRDS("ignore/prob_hosp/paton_inferred.rds") %>%
  select(pfpr, sma, distance, fever_tx, py, country)

# Load fit
parameters <- readRDS("ignore/prob_hosp/mcmc_fits/parameters.rds") %>%
  mutate(country = ifelse(country == "Congo Democratic Republic", "DRC", country)) 

# Prediction - DHS
predict <- parameters %>%
  slice_sample(n = 100) %>%
  select(sample, country, global_capacity, country_capacity, shift, pfpr_beta, distance_beta, tx_beta) %>%
  left_join(dhs_sma, by = "country") %>%
  mutate(predict = pmap_dbl(select(., -country, -sma_microscopy, -sample, - pfprg, -distanceg, -txg), gompertz)) %>%
  filter(!is.na(predict))


pd1 <- predict %>%
  mutate(ml = mean(sma_microscopy))

ppc1 <- ggplot(data = pd1, aes(x = predict)) +
  geom_histogram(fill = "darkblue", bins = 30) +
  geom_vline(aes(xintercept = ml), col = "red") +
  xlab("SMA prevalence") +
  theme_bw()

pd2 <- predict %>%
  group_by(country) %>%
  mutate(ml = mean(sma_microscopy))

ppc2 <- ggplot(data = pd2, aes(x = predict)) +
  geom_histogram(fill = "darkblue", bins = 20) +
  geom_vline(aes(xintercept = ml), col = "red") +
  xlab("SMA prevalence") +
  theme_bw() +
  facet_wrap(~country, ncol = 6, scales = "free") +
  theme(strip.background = element_rect(fill = NA))

pd3 <- predict %>%
  group_by(pfprg) %>%
  mutate(ml = mean(sma_microscopy))

ppc3 <- ggplot(data = pd3, aes(x = predict)) +
  geom_histogram(fill = "darkblue", bins = 20) +
  geom_vline(aes(xintercept = ml), col = "red") +
  xlab("SMA prevalence") +
  theme_bw() +
  facet_wrap(~pfprg, scales = "free") +
  theme(strip.background = element_rect(fill = NA))

pd4 <- predict %>%
  group_by(distanceg) %>%
  mutate(ml = mean(sma_microscopy))

ppc4 <- ggplot(data = pd4, aes(x = predict)) +
  geom_histogram(fill = "darkblue", bins = 20) +
  geom_vline(aes(xintercept = ml), col = "red") +
  xlab("SMA prevalence") +
  theme_bw() +
  facet_wrap(~distanceg, scales = "free") +
  theme(strip.background = element_rect(fill = NA))

pd5 <- predict %>%
  group_by(txg) %>%
  mutate(ml = mean(sma_microscopy))

ppc5 <- ggplot(data = pd5, aes(x = predict)) +
  geom_histogram(fill = "darkblue", bins = 20) +
  geom_vline(aes(xintercept = ml), col = "red") +
  xlab("SMA prevalence") +
  theme_bw() +
  facet_wrap(~txg, scales = "free") +
  theme(strip.background = element_rect(fill = NA))


# Prediction - Paton
predict2 <- parameters %>%
  slice_sample(n = 500)%>%
  filter(country %in% unique(paton$country)) %>%
  left_join(paton, by = "country") %>%
  mutate(sma_prevalence = pmap_dbl(select(., -country, -py, -prob_symptomatic, -dur, -prob_recognise, -hosp, -group_sd, -sma, -sample), gompertz),
         adjusted_sma_prevalence = sma_prev_age_standardise(sma_prevalence),
         symptomatic_sma_prevalence = adjusted_sma_prevalence * prob_symptomatic,
         community_symptomatic_sma_inc = inc1(prevalence = symptomatic_sma_prevalence, recovery_rate = 1 / dur),
         community_recognised_sma = community_symptomatic_sma_inc * prob_recognise,
         predict = community_recognised_sma * hosp) %>%
  filter(!is.na(predict))

pd6 <- predict2 %>%
  mutate(ml = mean(1000 * sma / py))

ppc6 <- ggplot(data = pd6, aes(x = predict)) +
  geom_histogram(fill = "darkblue", bins = 40) +
  geom_vline(aes(xintercept = ml), col = "red") +
  xlab("Hospitalised incidence per 1000") +
  theme_bw()

pd7 <- predict2 %>%
  group_by(country) %>%
  mutate(ml = mean(1000 * sma / py))

ppc7 <- ggplot(data = pd7, aes(x = predict)) +
  geom_histogram(fill = "darkblue", bins = 20) +
  geom_vline(aes(xintercept = ml), col = "red") +
  xlab("Hospitalised incidence per 1000") +
  theme_bw() +
  facet_wrap(~country, scales = "free") +
  theme(strip.background = element_rect(fill = NA))