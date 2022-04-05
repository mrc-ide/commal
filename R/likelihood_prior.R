r_loglike <- function(params, data, misc) {
  country_capacity <- params[grepl("ccc", names(params))]
  country_hosp <- params[grepl("hosp_", names(params))]
  block <- unlist(misc["block"])
  n_countries = misc$n_countries
  
  if(block %in% 1:n_countries){
    prob_dhs <- misc$model_function(pfpr = data$dhs[[block]]$pfpr,
                                    global_capacity = params["global_capacity"],
                                    country_capacity = country_capacity[block],
                                    pfpr_beta = params["pfpr_beta"],
                                    shift = params["shift"])
    
    loglike <- sum(dbinom(data$dhs[[block]]$symp_sma_microscopy, 1, prob_dhs, log = T)) +
      sum(dbinom(data$dhs[[block]]$chronic_anaemia, 1, params["chronic"], log = T))
  }
  
  if(block %in% (n_countries + 1):(n_countries + 3)){
    paton_block <- block - n_countries
    
    prob_paton <- misc$model_function(pfpr = data$paton[[paton_block]]$pfpr,
                                      global_capacity = params["global_capacity"],
                                      country_capacity = country_capacity[paton_block],
                                      pfpr_beta = params["pfpr_beta"],
                                      shift = params["shift"])
    
    hosp_inc <- cascade(symptomatic_sma_prevalence = prob_paton,
                        chronic = params["chronic"],
                        dur = params["dur"],
                        py = data$paton[[paton_block]]$py, 
                        prob_recognise = params["prob_recognise"],
                        hosp = country_hosp[paton_block],
                        distance_beta = params["distance_beta"],
                        distance = data$paton[[paton_block]]$distance)
    
    loglike <- sum(dnbinom(data$paton[[paton_block]]$sma, mu = hosp_inc, size = params["overdispersion"], log = T))
  }
  
  if(block == (n_countries + 4)){
    loglike <- sum(dnorm(country_capacity, 0, params["group_sd"], log = T))
  }
  
  # return
  return(loglike)
}

r_logprior <- function(params, misc){
  ret <- 
    sum(dnorm(params[c("global_capacity", "shift", "pfpr_beta")], 0, 10, log = TRUE)) +
    # Probably a minimum bound: Mean 4.61 of: from Mousa (2020) supplement data S1: filter(SMA = 1, age between 3 months and 9 years).
    sum(dgamma2(params["dur"], mean = 14, var = 150, log = TRUE)) +
    sum(dunif(params["chronic"], 0, 1, log = TRUE)) +
    # Prior from Shellenberg (2003) Figure 2 (see paragraph text)
    sum(dbeta(params["prob_recognise"], 14, 43 - 14, log = TRUE)) +
    sum(dlnorm(params["overdispersion"], 0, 5, log = TRUE)) +
    #sum(dlnorm(params["group_sd"], 0, 5, log = TRUE)) +
    sum(dunif(params["group_sd"], 0, 10000, log = TRUE)) +
    sum(dnorm(params[grepl("ccc", names(params))], 0, 10, log = TRUE)) +
    # Weakly informative prior on prop hospitalisation (mean set to match mean across all countries of mu_DA from Camponovo)
    #sum(dnorm(params[grepl("hosp_", names(params))], -0.35, 1, log = TRUE)) +
    sum(dnorm(params[grepl("hosp_", names(params))], 0, 1.4, log = TRUE)) +
    sum(dnorm(params["distance_beta"], 0, 10, log = TRUE))
  return(ret)
}
