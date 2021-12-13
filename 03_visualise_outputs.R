library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

#mcmc <- readRDS("ignore/prob_hosp/mcmc_fits/mcmc.rds")

# Extract estimate of best par
best_par <- mcmc$output %>%
  filter(phase == "sampling") %>%
  select(-c(chain, phase, iteration, logprior, loglikelihood)) %>%
  summarise_all(median) %>%  #slice_max(order_by = (logprior + loglikelihood), n = 1) %>%
  unlist()
samples <- sample_chains(mcmc, 500)

### Check asymptotes and country-level random effects ##########################
country_capacity <- best_par[grepl("ccc", names(best_par))] %>%
  t() %>%
  as_tibble() %>%
  pivot_longer(-c(), names_to = "country", values_to = "country_capacity", names_prefix = "ccc_") %>%
  mutate(asymptote = rlogit(best_par["global_capacity"] + country_capacity))
country_capacity_plot <- ggplot(country_capacity, aes(y = country_capacity, x = country)) +
  geom_hline(yintercept = c(0, 1.96 * best_par["group_sd"], -1.96 * best_par["groud_sd"]), lty = c(1, 2, 2)) +
  geom_point() +
  ylab("Country random effect") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())
asymptote_plot <- ggplot(country_capacity, aes(y = asymptote, x = country)) +
  geom_point() +
  ylab("Country asymptote (distance = 0)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())
################################################################################

### Compare pfpr fit to DHS data ###############################################
dhs_summary_pfpr <- dhs_sma %>%
  bind_rows() %>%
  group_by(country) %>%
  mutate(pfpr_g = cut_interval(pfpr, 7)) %>%
  group_by(country, pfpr_g) %>%
  summarise(pfpr = mean(pfpr),
            smal = binom::binom.exact(sum(sma), n())$lower,
            smau = binom::binom.exact(sum(sma), n())$upper,
            sma_prevalence = mean(sma),
            distance = mean(distance))

dhs_predict_pfpr <- dhs_summary_pfpr %>%
  group_by(country) %>%
  summarise(distance = mean(distance)) %>%
  left_join(data.frame(pfpr = seq(0, 0.7, 0.01)), by = character()) %>%
  left_join(country_capacity, by = "country") %>%
  mutate(predicted_sma_prevalence = misc$model_function(pfpr = pfpr,
                                                        distance = distance, 
                                                        global_capacity = best_par["global_capacity"],
                                                        country_capacity = country_capacity,
                                                        distance_beta = best_par["distance_beta"],
                                                        pfpr_beta = best_par["pfpr_beta"],
                                                        shift = best_par["shift"]))
pfpr_dhs_plot <- ggplot() + 
  geom_point(data = dhs_summary_pfpr, aes(x = pfpr, y = sma_prevalence)) +
  geom_linerange(data = dhs_summary_pfpr, aes(x = pfpr, y = sma_prevalence, ymin = smal, ymax = smau)) + 
  geom_line(data = dhs_predict_pfpr, aes(x = pfpr, y = predicted_sma_prevalence)) +
  xlab("PfPr") +
  ylab("SMA prevalence") +
  theme_bw() +
  ggtitle("Model fit to DHS data") +
  theme(strip.background = element_rect(fill = NA)) +
  facet_wrap(~ country, scales = "free_y")

pfpr_sample_pd <- samples %>%
  select(sample, global_capacity, shift, pfpr_beta) %>%
  left_join(data.frame(pfpr = seq(0, 0.7, 0.01)), by = character()) %>%
  mutate(predicted_sma_prevalence = misc$model_function(pfpr = pfpr,
                                                        distance = 1, 
                                                        global_capacity = global_capacity,
                                                        country_capacity = 0,
                                                        distance_beta = 0,
                                                        pfpr_beta = pfpr_beta,
                                                        shift = shift))
pfpf_best_pd <- best_par %>%
  t() %>%
  as_tibble() %>%
  select(global_capacity, shift, pfpr_beta) %>%
  left_join(data.frame(pfpr = seq(0, 0.7, 0.01)), by = character()) %>%
  mutate(predicted_sma_prevalence = misc$model_function(pfpr = pfpr,
                                                        distance = 1, 
                                                        global_capacity = global_capacity,
                                                        country_capacity = 0,
                                                        distance_beta = 0,
                                                        pfpr_beta = pfpr_beta,
                                                        shift = shift))

pfpr_sample_plot <- ggplot() +
  geom_line(data = pfpr_sample_pd, aes(x = pfpr, y = predicted_sma_prevalence, group = sample), alpha = 0.1) +
  geom_line(data = pfpf_best_pd, aes(x = pfpr, y = predicted_sma_prevalence), col = "dodgerblue", size = 1) +
  xlab("PfPr") +
  ylab("SMA prevalence") +
  theme_bw()
################################################################################

## Compare distance fit to DHS data #################################################
dhs_summary_distance <- dhs_sma %>%
  bind_rows() %>%
  group_by(country) %>%
  mutate(distance_g = cut_interval(distance, 7)) %>%
  group_by(country, distance_g) %>%
  summarise(pfpr = mean(pfpr),
            smal = binom::binom.exact(sum(sma), n())$lower,
            smau = binom::binom.exact(sum(sma), n())$upper,
            sma_prevalence = mean(sma),
            distance = mean(distance))

max_distance <- dhs_summary_distance %>%
  group_by(country) %>%
  slice_max(distance) %>%
  select(country, distance) %>%
  rename(max_distance = distance)

dhs_predict_distance <- dhs_summary_pfpr %>%
  group_by(country) %>%
  summarise(pfpr = mean(pfpr)) %>%
  left_join(data.frame(distance = 0:500), by = character()) %>%
  left_join(max_distance, by = "country") %>%
  filter(distance < max_distance) %>%
  left_join(country_capacity, by = "country") %>%
  mutate(predicted_sma_prevalence = misc$model_function(pfpr = pfpr,
                                                        distance = distance, 
                                                        global_capacity = best_par["global_capacity"],
                                                        country_capacity = country_capacity,
                                                        distance_beta = best_par["distance_beta"],
                                                        pfpr_beta = best_par["pfpr_beta"],
                                                        shift = best_par["shift"]))
distance_dhs_plot <- ggplot() + 
  geom_point(data = dhs_summary_distance, aes(x = distance, y = sma_prevalence)) +
  geom_linerange(data = dhs_summary_distance, aes(x = distance, y = sma_prevalence, ymin = smal, ymax = smau)) + 
  geom_line(data = dhs_predict_distance, aes(x = distance, y = predicted_sma_prevalence)) +
  xlab("Distance") +
  ylab("SMA prevalence") +
  theme_bw() +
  ggtitle("Model fit to DHS data") +
  theme(strip.background = element_rect(fill = NA)) +
  facet_wrap(~ country, scales = "free")


distance_sample_pd <- samples %>%
  select(sample, global_capacity, shift, distance_beta, pfpr_beta) %>%
  left_join(data.frame(distance = 0:200), by = character()) %>%
  mutate(predicted_sma_prevalence = misc$model_function(pfpr = 0.2,
                                                        distance = distance, 
                                                        global_capacity = global_capacity,
                                                        country_capacity = 0,
                                                        distance_beta = distance_beta,
                                                        pfpr_beta = pfpr_beta,
                                                        shift = shift))
distance_best_pd <- best_par %>%
  t() %>%
  as_tibble() %>%
  select(global_capacity, shift, distance_beta, pfpr_beta) %>%
  left_join(data.frame(distance = 0:200), by = character()) %>%
  mutate(predicted_sma_prevalence = misc$model_function(pfpr = 0.2,
                                                        distance = distance, 
                                                        global_capacity = global_capacity,
                                                        country_capacity = 0,
                                                        distance_beta = distance_beta,
                                                        pfpr_beta = pfpr_beta,
                                                        shift = shift))

distance_sample_plot <- ggplot() +
  geom_line(data = distance_sample_pd, aes(x = distance, y = predicted_sma_prevalence, group = sample), alpha = 0.1) +
  geom_line(data = distance_best_pd, aes(x = distance, y = predicted_sma_prevalence), col = "dodgerblue", size = 1) +
  xlab("Distance") +
  ylab("SMA prevalence") +
  theme_bw()
################################################################################

### Compare output to Paton et al hospital data ################################
paton_hosp <- best_par %>%
  t() %>%
  as_tibble() %>%
  select(contains("hosp")) %>%
  pivot_longer(-c(), names_to = "country", values_to = "hosp", names_prefix = "hosp_") %>%
  rename(hospor = hosp)

paton_average <- paton %>%
  bind_rows() %>%
  group_by(country) %>%
  summarise(
    distance = mean(distance),
    act = mean(act)
  ) %>%
  left_join(paton_hosp, by = "country") %>%
  mutate(cfr = best_par["cfr"],
         dur_recover = best_par["dur_recover"],
         dur_tx  =  best_par["dur_tx"],
         dur_die = best_par["dur_die"],
         p_recover = (1 - act) * (1 -cfr),
         p_tx =   act,
         p_die = (1 - act) * cfr) %>%
  rowwise() %>%
  mutate(rate = 1 / weighted.mean(c(dur_recover, dur_tx, dur_die), c(p_recover, p_tx, p_die)))

paton_fit <- data.frame(
  country = c("Kenya", "Tanzania", "Uganda")) %>%
  left_join(data.frame(pfpr = seq(0, 0.7, 0.01)), by = character()) %>%
  left_join(country_capacity, by = "country") %>%
  left_join(paton_average, by = "country") %>%
  mutate(predicted_sma_prevalence = misc$model_function(pfpr = pfpr,
                                                        distance = distance, 
                                                        global_capacity = best_par["global_capacity"],
                                                        country_capacity = country_capacity,
                                                        distance_beta = best_par["distance_beta"],
                                                        pfpr_beta = best_par["pfpr_beta"],
                                                        shift = best_par["shift"])) %>%
  mutate(predicted_sma_prevalence = sma_prev_age_standardise(predicted_sma_prevalence),
         community = 1000 *  365 * inc1(predicted_sma_prevalence, rate),
         hospital = community / hospor,
         paton_average = paton_sma(pfpr))

hosp_community_pd <- paton_fit %>%
  select(country, pfpr, community, hospital) %>%
  pivot_longer(-c(country, pfpr), names_to = "location", values_to = "inc")

paton_fit_plot <- ggplot() +
  geom_point(data = bind_rows(paton), aes(x = pfpr, y = 1000 * (sma / py))) +
  geom_line(data = paton_fit, aes(x = pfpr, y = paton_average), lty = 2) +
  geom_line(data = paton_fit, aes(x = pfpr, y = hospital)) +
  xlab("PfPr") +
  ylab("Annual hospitalised incidence per 1000 children") +
  theme_bw() +
  ggtitle("Model fit to Paton data") +
  facet_wrap(~ country) +
  theme(strip.background = element_rect(fill = NA))

hosp_community_inc_plot <- ggplot(hosp_community_pd, aes(x = pfpr, y = inc, col = location)) +
  geom_line() +
  xlab("PfPr") +
  ylab("Annual incidence per 1000 children") +
  theme_bw() +
  ggtitle("Model fit to Paton data") +
  facet_wrap(~ country) +
  theme(strip.background = element_rect(fill = NA))
################################################################################

### Probability of hospitalisation #############################################
prob_hosp <- mcmc$output %>%
  filter(phase == "sampling") %>%
  select(contains("hosp")) %>%
  tidyr::pivot_longer(-c(), names_to = "country", values_to = "hosp_or", names_prefix = "hosp_") %>%
  mutate(prop_hosp = 1 - p(hosp_or))

prob_hosp_plot <- ggplot(prob_hosp, aes(x = prop_hosp)) +
  geom_histogram(binwidth = 0.01) + 
  facet_wrap( ~ country) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA))

prob_hosp %>%
  group_by(country) %>%
  summarise(prop_hosp_l = quantile(prop_hosp, 0.025),
            prop_hosp = quantile(prop_hosp, 0.5),
            prop_hosp_u = quantile(prop_hosp, 0.975))
  
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