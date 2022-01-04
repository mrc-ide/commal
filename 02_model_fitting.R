# Load packages
library(drjacoby)
library(dplyr)
library(ggplot2)
library(patchwork)

# Load functions
source("R/prevalence_incidence.R")
source("R/gamma.R")
source("R/odds_probability.R")
source("R/paton_fits.R")
source("R/model_functions.R")
source("R/likelihood_prior.R")

# Load data
paton <- readRDS("ignore/prob_hosp/paton_inferred.rds")%>%
  select(pfpr, sma, distance, act, py, country) %>%
  split(.$country)
paton_countries <-  as.character(sapply(paton, function(x)x$country[1]))

dhs_sma <- readRDS("ignore/prob_hosp/dhs_sma.rds") %>%
  select(pfpr, sma_microscopy, distance, country) %>%
  split(.$country)
cn <- as.character(sapply(dhs_sma, function(x)x$country[1]))
names(dhs_sma) <- cn
reorder <- c(paton_countries, setdiff(cn, paton_countries))
dhs_sma <- dhs_sma[reorder]

country_names <- names(dhs_sma)
n_countries <- length(dhs_sma)

# Data input list for MCMC
data_list <- list(
  dhs = lapply(dhs_sma, function(x) x[,c("pfpr", "distance", "sma_microscopy")]),
  paton = lapply(paton, function(x) x[,c("pfpr", "distance", "act", "py", "sma")])
)

# Helper functions for MCMC
misc <- list(
  model_function = gompertz,
  sma_prev_age_standardise = sma_prev_age_standardise,
  mean_duration = mean_duration
)

# Define parameters
global_params <- define_params(name = "global_capacity", min = -Inf, max = Inf, init = -6, block = 1:(n_countries+3), 
                               name = "shift", min = 0, max = 50, init = 5, block = 1:(n_countries+3), 
                               name = "pfpr_beta", min = 0, max = 50, init = 10, block = 1:(n_countries+3), 
                               name = "distance_beta", min = -Inf, max = Inf, init = 0, block = 1:(n_countries+3),
                               name = "dur", min = 0, max = Inf, init = 365, block = (n_countries+1):(n_countries+3),
                               name = "group_sd", min = 0, max = Inf, init = 1, block = (n_countries+4))
country_params <- data.frame(
  name = paste0("ccc_", country_names),
  min = -10,
  max = 10
)
country_params$init <- as.list(rep(0, n_countries))
country_params$block <- lapply(1:n_countries, function(x){
  c(x, (n_countries+4))
})
country_params$block[[1]] <- c(1, (n_countries+1), (n_countries+4))
country_params$block[[2]] <- c(2, (n_countries+2), (n_countries+4))
country_params$block[[3]] <- c(3, (n_countries+3), (n_countries+4))

hosp_params <- data.frame(
  name = paste0("hosp_", paton_countries),
  min = 0,
  max = 1000
)
hosp_params$init <- as.list(rep(1, 3))
hosp_params$block <- as.list((n_countries+1):(n_countries+3))
df_params <- bind_rows(global_params, hosp_params, country_params)

# Run MCMC
#cl <- parallel::makeCluster(4)
#parallel::clusterExport(cl, c("rlogit", "inc1", "dgamma2"))
mcmc <- run_mcmc(data = data_list,
                 df_params = df_params,
                 loglike = r_loglike,
                 logprior = r_logprior,
                 misc = misc,
                 burnin = 5e3,
                 samples = 5e3,
                 rungs = 1,
                 chains = 1)
#parallel::stopCluster(cl)
#saveRDS(mcmc, "ignore/prob_hosp/mcmc_fits/mcmc.rds")


