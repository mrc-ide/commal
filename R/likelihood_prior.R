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
    
    loglike <- sum(dbinom(data$dhs[[block]]$sma, 1, prob_dhs, log = T))
  }
  
  if(block %in% (n_countries + 1):(n_countries + 3)){
    chronic <- params[grepl("chronic", names(params))] 
    
    paton_block <- block - n_countries
    
    prob_paton <- misc$model_function(pfpr = data$paton[[paton_block]]$pfpr,
                                      global_capacity = params["global_capacity"],
                                      country_capacity = country_capacity[paton_block],
                                      pfpr_beta = params["pfpr_beta"],
                                      shift = params["shift"])
    
    sa_input <- 0
    if(misc$adjust_sa){
      sa_input <- chronic[paton_block]
    }
    
    hosp_inc <- cascade(sma_prevalence = prob_paton,
                        chronic = sa_input,
                        pfpr = data$paton[[paton_block]]$pfpr,
                        dur = params["dur"],
                        py = data$paton[[paton_block]]$py, 
                        hosp = country_hosp[paton_block],
                        distance_beta = params["distance_beta"],
                        distance = data$paton[[paton_block]]$distance)
    
    loglike <- sum(dnbinom(data$paton[[paton_block]]$sma, mu = hosp_inc, size = params["overdispersion"], log = T)) +
      sum(dbinom(data$dhs[[paton_block]]$sa[data$dhs[[paton_block]]$diagnostic == "negative"], 1, chronic[paton_block], log = T))
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
    # Setting duration with a mean of 14 days and upper 95% quantile of 60 days
    # Probably a minimum bound: Mean 4.61 of: from Mousa (2020) supplement data S1: filter(SMA = 1, age between 3 months and 9 years).
    sum(dlnorm(params["dur"], log(14), 0.7425, log = TRUE)) +
    #sum(dgamma(params["dur"], shape = 0.71, rate = 0.0507, log = TRUE)) +
    #sum(dgamma(params["dur"], shape = 1.546, rate = 0.335, log = TRUE)) +
    sum(dunif(params[grepl("chronic", names(params))], 0, 1, log = TRUE)) +
    sum(dlnorm(params["overdispersion"], 0, 5, log = TRUE)) +
    sum(dunif(params["group_sd"], 0, 10000, log = TRUE)) +
    sum(dnorm(params[grepl("ccc", names(params))], 0, 10, log = TRUE)) +
    sum(dnorm(params[grepl("hosp_", names(params))], 0, 1.5, log = TRUE)) +
    sum(dnorm(params["distance_beta"], 0, 10, log = TRUE))
  return(ret)
}

r_logprior_dur <- function(params, misc){
  ret <- 
    sum(dnorm(params[c("global_capacity", "shift", "pfpr_beta")], 0, 10, log = TRUE)) +
    # Setting duration with a mean of 14 days and upper 95% quantile of ~365 days
    sum(dlnorm(params["dur"], log(14), 1.66, log = TRUE)) +
    sum(dunif(params[grepl("chronic", names(params))], 0, 1, log = TRUE)) +
    sum(dlnorm(params["overdispersion"], 0, 5, log = TRUE)) +
    sum(dunif(params["group_sd"], 0, 10000, log = TRUE)) +
    sum(dnorm(params[grepl("ccc", names(params))], 0, 10, log = TRUE)) +
    sum(dnorm(params[grepl("hosp_", names(params))], 0, 1.5, log = TRUE)) +
    sum(dnorm(params["distance_beta"], 0, 10, log = TRUE))
  return(ret)
}
