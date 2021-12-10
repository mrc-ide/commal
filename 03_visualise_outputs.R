library(ggplot2)
library(patchwork)

mcmc <- readRDS("ignore/prob_hosp/mcmc_fits/mcmc.rds")

# Extract estimate of best par
best_par <- mcmc$output %>%
  filter(phase == "sampling") %>%
  select(-c(chain, phase, iteration, logprior, loglikelihood)) %>%
  summarise_all(median) %>%
  unlist()

## Compare GLM fit to DHS data #################################################
ccc <- best_par[grepl("ccc", names(best_par))]
plot(rlogit(best_par["a"] + ccc), ylab = "Asymptotes")
plot(ccc, ylim = c(-2, 2))
abline(h = 0)
abline(h = 1.96 * best_par["group_sd"], lty = 2)
abline(h = -1.96 * best_par["group_sd"], lty = 2)

dhs_pd <- list()
dhs_predict_pd <- list()
for(i in 1:length(dhs_sma)){
  dhs_pd[[i]] <- dhs_sma[[i]] %>%
    mutate(pfpr_g = cut_number(pfpr, 5)) %>%
    group_by(country, pfpr_g) %>%
    summarise(pfpr = mean(pfpr),
              smal = binom::binom.exact(sum(sma), n_countries())$lower,
              smau = binom::binom.exact(sum(sma), n_countries())$upper,
              sma_prevalence = mean(sma),
              distance = mean(distance))
  
  dhs_predict_pd[[i]] <- data.frame(
    pfpr = seq(0, 0.6, 0.01),
    country = country_names[i],
    distance = mean(dhs_pd[[i]]$distance)) %>%
    mutate(predicted_sma_prevalence = misc$model_function(pfpr, distance, best_par["a"], best_par["b"], best_par["c"], ccc[i], best_par["e"]))
}
dhs_pd <- bind_rows(dhs_pd)
dhs_predict_pd <- bind_rows(dhs_predict_pd)

ggplot() + 
  geom_point(data = dhs_pd, aes(x = pfpr, y = sma_prevalence)) +
  geom_linerange(data = dhs_pd, aes(x = pfpr, y = sma_prevalence, ymin = smal, ymax = smau)) + 
  geom_line(data = dhs_predict_pd, aes(x = pfpr, y = predicted_sma_prevalence)) +
  xlab("PfPr") +
  ylab("SMA prevalence") +
  theme_bw() +
  ggtitle("Model fit to DHS data") +
  theme(strip.background = element_rect(fill = NA)) +
  facet_wrap(~country)
################################################################################

## Compare GLM fit to disatne data #################################################

dhs_pd <- list()
dhs_predict_pd <- list()
for(i in 1:length(dhs_sma)){
  dhs_pd[[i]] <- dhs_sma[[i]] %>%
    mutate(distance_g = cut_number(distance, 10)) %>%
    group_by(country, distance_g) %>%
    summarise(pfpr = mean(pfpr),
              smal = binom::binom.exact(sum(sma), n_countries())$lower,
              smau = binom::binom.exact(sum(sma), n_countries())$upper,
              sma_prevalence = mean(sma),
              distance = mean(distance))
  
  dhs_predict_pd[[i]] <- data.frame(
    distance = seq(0, max(dhs_pd[[i]]$distance), length.out = 100),
    country = country_names[i],
    pfpr = mean(dhs_pd[[i]]$pfpr)) %>%
    mutate(predicted_sma_prevalence = misc$model_function(pfpr, distance, best_par["a"], best_par["b"], best_par["c"], ccc[i], best_par["e"]))
}
dhs_pd <- bind_rows(dhs_pd)
dhs_predict_pd <- bind_rows(dhs_predict_pd)

ggplot() + 
  geom_point(data = dhs_pd, aes(x = distance, y = sma_prevalence)) +
  geom_linerange(data = dhs_pd, aes(x = distance, y = sma_prevalence, ymin = smal, ymax = smau)) + 
  geom_line(data = dhs_predict_pd, aes(x = distance, y = predicted_sma_prevalence)) +
  xlab("Distance") +
  ylab("SMA prevalence") +
  theme_bw() +
  ggtitle("Model fit to DHS data") +
  theme(strip.background = element_rect(fill = NA)) +
  facet_wrap(~country, scales = "free_x")
################################################################################