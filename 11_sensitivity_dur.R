### Model fitting ##############################################################

# Load packages
library(drjacoby)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)

# Load functions
source("R/prevalence_incidence.R")
source("R/gamma.R")
source("R/odds_probability.R")
source("R/paton_fits.R")
source("R/model_functions.R")
source("R/likelihood_prior.R")
source("R/cascade.R")

# Load data
paton <- readRDS("ignore/prob_hosp/paton_inferred.rds") %>%
  rename(sma = sma_n_diamond) %>%
  select(pfpr, sma, distance, py, country) %>%
  split(.$country)
paton_countries <-  as.character(sapply(paton, function(x)x$country[1]))

dhs_masa <- readRDS("ignore/prob_hosp/dhs_masa.rds") %>%
  mutate(diagnostic = microscopy) %>%
  select(pfpr, masa, sa, diagnostic, country) %>%
  split(.$country)
cn <- as.character(sapply(dhs_masa, function(x)x$country[1]))
names(dhs_masa) <- cn
reorder <- c(paton_countries, setdiff(cn, paton_countries))
dhs_masa <- dhs_masa[reorder]

country_names <- names(dhs_masa)
n_countries <- length(dhs_masa)

# Data input list for MCMC
data_list <- list(
  dhs = lapply(dhs_masa, function(x) x[,c("pfpr", "masa", "sa", "diagnostic")]),
  paton = lapply(paton, function(x) x[,c("pfpr", "distance", "py", "sma")])
)

# Helper functions for MCMC
misc <- list(
  model_function = gompertz,
  n_countries = n_countries,
  adjust_sa = TRUE
)

################################################################################
### Model 1 ####################################################################
################################################################################

# Define parameters
global_params <- define_params(name = "global_capacity", min = -Inf, max = Inf, init = c(-7, -6, -5, -4), block = 1:(n_countries+3), 
                               name = "shift", min = 0, max = Inf, init = c(1, 1.3, 1.2, 1.3), block = 1:(n_countries+3), 
                               name = "pfpr_beta", min = -Inf, max = Inf, init = 6:9, block = 1:(n_countries+3), 
                               name = "dur", min = 0, max = Inf, init = c(2, 4, 6, 8), block = (n_countries+1):(n_countries+3),
                               name = "distance_beta", min = -Inf, max = Inf, init = c(-0.1, -0.05, 0.05, 0.1), block = (n_countries+1):(n_countries+3), 
                               name = "overdispersion", min = 0, max = Inf, init = 1:4, block = (n_countries+1):(n_countries+3),
                               name = "group_sd", min = 0, max = Inf, init = c(0.5, 0.75, 1, 1.25), block = (n_countries+4))
country_params <- data.frame(
  name = paste0("ccc_", country_names),
  min = -Inf,
  max = Inf
)
country_params$init <- lapply(1:n_countries, function(x){
  c(-0.1, -0.05, 0.05, 0.1)
})
country_params$block <- lapply(1:n_countries, function(x){
  c(x, (n_countries+4))
})
country_params$block[[1]] <- c(1, (n_countries+1), (n_countries+4))
country_params$block[[2]] <- c(2, (n_countries+2), (n_countries+4))
country_params$block[[3]] <- c(3, (n_countries+3), (n_countries+4))

chronic_params <- data.frame(
  name = paste0("chronic_", paton_countries),
  min = 0,
  max = 1
)
chronic_params$init <- lapply(1:3, function(x){
  c(0.0001, 0.0002, 0.0003, 0.0004)
})
chronic_params$block <- as.list((n_countries+1):(n_countries+3))

hosp_params <- data.frame(
  name = paste0("hosp_", paton_countries),
  min = -Inf,
  max = Inf
)
hosp_params$init <- lapply(1:3, function(x){
  c(-0.02, -0.01, 0.01, 0.02)
})
hosp_params$block <- as.list((n_countries+1):(n_countries+3))

df_params <- bind_rows(global_params, chronic_params, hosp_params, country_params)

# Run MCMC
cl <- parallel::makeCluster(4)
parallel::clusterExport(cl, c("rlogit", "prev_to_inc", "dgamma2",
                              "cascade", "%>%", "malaria_attributable",
                              "sma_prev_age_standardise", "hospitalised"))
mcmc <- run_mcmc(data = data_list,
                 df_params = df_params,
                 loglike = r_loglike,
                 logprior = r_logprior_dur,
                 misc = misc,
                 burnin = 50000,
                 samples = 50000,
                 rungs = 1,
                 chains = 4,
                 cluster = cl)
parallel::stopCluster(cl)
saveRDS(mcmc, "ignore/prob_hosp/mcmc_fits/mcmc_sensitivity_dur.rds")
#plot_par(mcmc)

################################################################################
################################################################################
################################################################################

################################################################################
### Wrangle parameters #########################################################
################################################################################
samples <- mcmc$output %>%
  filter(phase == "sampling") %>%
  select(-c(chain, phase, iteration, logprior, loglikelihood)) %>%
  mutate(sample = 1:n())

global_parameters <- samples %>%
  select(-contains("ccc_"), -contains("hosp"), -contains("chronic"))

country_parameters <- samples %>%
  select(contains("ccc_"), sample) %>%
  pivot_longer(-sample, names_to = "country", values_to = "country_capacity", names_prefix = "ccc_") 

chronic_parameters <- samples %>%
  select(contains("chronic_"), sample) %>%
  pivot_longer(-sample, names_to = "country", values_to = "chronic", names_prefix = "chronic_") 

hospital_parameters <- samples %>%
  select(contains("hosp_"), sample) %>%
  pivot_longer(-sample, names_to = "country", values_to = "hosp", names_prefix = "hosp_") 

parameters <- global_parameters %>%
  left_join(country_parameters, by = "sample") %>%
  left_join(chronic_parameters, by = c("country", "sample")) %>%
  left_join(hospital_parameters, by = c("country", "sample")) %>%
  select(country, sample, everything())

saveRDS(parameters, "ignore/prob_hosp/mcmc_fits/parameters_sensitivity_dur.rds")
################################################################################
################################################################################
################################################################################