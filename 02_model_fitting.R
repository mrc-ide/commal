# Load packages
library(drjacoby)
library(dplyr)
library(patchwork)


# Load functions
source("R/prevalence_incidence.R")
source("R/gamma.R")
source("R/odds_probability.R")
source("R/paton_fits.R")
source("R/model_functions.R")
source("R/likelihood_prior.R")

# Load data
paton <- readRDS("ignore/prob_hosp/paton_inferred.rds")
dhs_sma <- readRDS("ignore/prob_hosp/dhs_sma.rds") %>%
  select(pfpr, sma, distance, country) %>%
  split(.$country)

country_names <- as.character(sapply(dhs_sma, function(x)x$country[1]))
n_countries <- length(dhs_sma)

# Data input list for MCMC
data_list <- lapply(dhs_sma, function(x) x[,c("pfpr", "distance", "sma")])

# Helper functions for MCMC
misc <- list(
  model_function = gompertz
)

# Define parameters
global_params <- define_params(name = "a", min = -Inf, max = Inf, init = -6, block = 1:n_countries, 
                           name = "b", min = 0, max = 50, init = 5, block = 1:n_countries, 
                           name = "c", min = 0, max = 50, init = 10, block = 1:n_countries, 
                           name = "e", min = -Inf, max = Inf, init = 0, block = 1:n_countries, 
                           name = "group_sd", min = 0, max = Inf, init = 1, block = n_countries+1)
country_params <- data.frame(
  name = paste0("ccc_", country_names),
  min = -10,
  max = 10
)
country_params$init <- as.list(rep(0, n_countries))
country_params$block <- lapply(1:n_countries, function(x){
  c(x, 23)
})
df_params <- bind_rows(global_params, country_params)

# Run MCMC
#cl <- parallel::makeCluster(2)
#parallel::clusterExport(cl, c("rlogit"))
mcmc <- run_mcmc(data = data_list,
                 df_params = df_params,
                 loglike = r_loglike,
                 logprior = r_logprior,
                 misc = misc,
                 burnin = 1e3,
                 samples = 1e4,
                 rungs = 1,
                 chains = 1)
#parallel::stopCluster(cl)
saveRDS(mcmc, "ignore/prob_hosp/mcmc_fits/mcmc.rds")

# Global parameters
pp <- plot_par(mcmc, c("a", "b", "c", "e", "group_sd"), display = FALSE)
trace_plots <- patchwork::wrap_plots(lapply(pp, function(x) x$trace + theme(legend.position = "none")), ncol = 5)
hist_plots <- patchwork::wrap_plots(lapply(pp, function(x) x$hist), ncol = 5)
acf_plots <- patchwork::wrap_plots(lapply(pp, function(x) x$acf), ncol = 5)
parameter_plots <- trace_plots / hist_plots / acf_plots

# Correlations between global pars
cor <- apply(combn(c("a", "b", "c", "e"), 2), 2, function(x){
  plot_cor(mcmc, x[1], x[2])
})
correlation_plots <- patchwork::wrap_plots(cor)
