r_loglike <- function(params, data, misc) {
  country_capacity <- params[grepl("ccc", names(params))]
  country_hosp <-params[grepl("hosp_", names(params))]
  block <- unlist(misc["block"])
  n_countries = misc$n_countries
  
  if(block %in% 1:n_countries){
    prob_dhs <- misc$model_function(pfpr = data$dhs[[block]]$pfpr,
                                    global_capacity = params["global_capacity"],
                                    country_capacity = country_capacity[block],
                                    pfpr_beta = params["pfpr_beta"],
                                    shift = params["shift"])
    loglike <- sum(dbinom(data$dhs[[block]]$sma_microscopy, 1, prob_dhs, log = T))
  }
  
  if(block %in% (n_countries + 1):(n_countries + 3)){
    paton_block <- block - n_countries
    prob_paton <- misc$model_function(pfpr = data$paton[[paton_block]]$pfpr,
                                      global_capacity = params["global_capacity"],
                                      country_capacity = country_capacity[paton_block],
                                      pfpr_beta = params["pfpr_beta"],
                                      shift = params["shift"])
      
    prob_paton <- misc$sma_prev_age_standardise(prob_paton)
    prob_paton <- prob_paton * params["prob_symptomatic"]
    inc <- inc1(prevalence = prob_paton,
                recovery_rate = 1 / params["dur"],
                py = data$paton[[paton_block]]$py)
    recognised_inc <- inc * params["prob_recognise"]
    hosp_inc <- recognised_inc * country_hosp[paton_block]
    loglike <- sum(dnbinom(data$paton[[paton_block]]$sma, mu = hosp_inc, size = 0.95, log = T))
  }
  
  if(block == (n_countries + 4)){
    loglike <- sum(dnorm(country_capacity, 0, params["group_sd"], log = T))
  }
  
  # return
  return(loglike)
}

r_logprior <- function(params, misc){
  ret <- 
    sum(dlnorm(params["group_sd"], 0, 5, log = TRUE)) +
    sum(dnorm(params[grepl("ccc", names(params))], 0, 10, log = TRUE)) +
    # Prior from Shellenberg (2003) Table 2
    sum(dbeta(params["prob_symptomatic"], 29, 64 - 29, log = TRUE)) +
    # Prior from Shellenberg (2003) Figure 2
    sum(dbeta(params["prob_recognise"], 14, 43 - 14, log = TRUE)) +
    sum(dlnorm(params["dur"], 0, 1, log = TRUE))
  return(ret)
}
