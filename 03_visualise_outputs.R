library(ggplot2)
library(patchwork)

#mcmc <- readRDS("ignore/prob_hosp/mcmc_fits/mcmc.rds")

# Extract estimate of best par
best_par <- mcmc$output %>%
  filter(phase == "sampling") %>%
  slice_max(order_by = (logprior + loglikelihood), n = 1) %>%
  select(-c(chain, phase, iteration, logprior, loglikelihood)) %>%
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
              smal = binom::binom.exact(sum(sma), n())$lower,
              smau = binom::binom.exact(sum(sma), n())$upper,
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
              smal = binom::binom.exact(sum(sma), n())$lower,
              smau = binom::binom.exact(sum(sma), n())$upper,
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


### Compare output to Paton et al hospital data ################################
act_cov <- mean(paton$act)
cfr <- as.numeric(best_par["cfr"])
probs <- c(
  p_recover = (1 - act_cov) * (1 -cfr),
  p_tx =   act_cov,
  p_die = (1 - act_cov) *cfr
)
ccc_paton <- ccc[paste0("ccc_", c("Kenya", "Tanzania", "Uganda"))]
mean_dist_paton <- paton %>%
  group_by(country) %>%
  summarise(distance = mean(distance),
            act = mean(act)) %>%
  mutate(countryn = case_when(
    country == "Kenya" ~ 1,
    country == "Tanzania" ~ 2,
    country == "Uganda" ~ 3
  ))

paton_fit <- expand.grid(
  pfpr = seq(0, 0.6, 0.01),
  countryn = 1:3) %>%
  left_join(mean_dist_paton, by = "countryn") %>%
  mutate(
    prob = misc$model_function(pfpr, distance, best_par["a"], best_par["b"], best_par["c"], ccc_paton[countryn], best_par["e"]),
    prob = misc$sma_prev_age_standardise(prob),
    p_recover = (1 - act) * (1 -cfr),
    p_tx = act,
    p_die = (1 - act) * cfr)
paton_fit$community <- apply(paton_fit[,6:9], 1, function(x, dur_recover, dur_tx, dur_die){
  1000 * 365 * inc1(x["prob"], 1 / weighted.mean(c(dur_recover, dur_tx, dur_die),
                                                 c(x["p_recover"], x["p_tx"], x["p_die"])))
}, dur_recover = best_par["dur_recover"], dur_tx = best_par["dur_tx"], dur_die = best_par["dur_die"])
paton_fit <- paton_fit %>%
  mutate(
    hosp = community / best_par[c("hosp_Kenya",  "hosp_Tanzania", "hosp_Uganda" )][countryn],
    paton_fit = paton_sma(pfpr)
  )

dp2 <- ggplot() +
  geom_point(data = paton, aes(x = pfpr, y = 1000 * (sma / py))) +
  geom_line(data = paton_fit, aes(x = pfpr, y = hosp), col = "black") +
  geom_line(data = paton_fit, aes(x = pfpr, y = paton_fit), col = "black" , lty = 2) +
  xlab("PfPr") +
  ylab("Annual hospitalised incidence per 1000 children") +
  theme_bw() +
  ggtitle("Model fit to Paton data") +
  facet_wrap(~ country) +
  theme(strip.background = element_rect(fill = NA))

dp3 <- ggplot() +
  geom_line(data = paton_fit, aes(x = pfpr, y = hosp), col = "red") +
  geom_line(data = paton_fit, aes(x = pfpr, y = community), col = "blue") +
  theme_bw() +
  ylab("Annual incidence per 1000 children") +
  ggtitle("Hospitalised vs community\nincidence") +
  facet_wrap(~country) +
  theme(strip.background = element_rect(fill = NA))
################################################################################