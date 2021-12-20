r_loglike <- function(params, data, misc) {
  country_capacity <- params[grepl("ccc", names(params))]
  country_hosp <-params[grepl("hosp_", names(params))]
  block <- unlist(misc["block"])
  
  if(block %in% 1:18){
    prob_dhs <- misc$model_function(pfpr = data$dhs[[block]]$pfpr,
                                    distance = data$dhs[[block]]$distance,
                                    global_capacity = params["global_capacity"],
                                    country_capacity = country_capacity[block],
                                    distance_beta = params["distance_beta"],
                                    pfpr_beta = params["pfpr_beta"],
                                    shift = params["shift"])
    loglike <- sum(dbinom(data$dhs[[block]]$sma_microscopy, 1, prob_dhs, log = T))
  }
  
  if(block %in% 19:21){
    paton_block <- block - 18
    prob_paton <- misc$model_function(pfpr = data$paton[[paton_block]]$pfpr,
                                      distance = data$paton[[paton_block]]$distance,
                                      global_capacity = params["global_capacity"],
                                      country_capacity = country_capacity[paton_block],
                                      distance_beta = params["distance_beta"],
                                      pfpr_beta = params["pfpr_beta"],
                                      shift = params["shift"])
    prob_paton <- misc$sma_prev_age_standardise(prob_paton)
    duration <- misc$mean_duration(dur_die = params["dur_die"],
                                   dur_recover = params["dur_recover"],
                                   cfr = params["cfr"])
    inc <- inc1(prevalence = prob_paton,
                recovery_rate = 1 / duration,
                py = data$paton[[paton_block]]$py)
    hosp_inc <- inc / country_hosp[paton_block]
    loglike <- sum(dpois(data$paton[[paton_block]]$sma, hosp_inc, log = T))
  }
  
  if(block == 22){
    loglike <- sum(dnorm(country_capacity, 0, params["group_sd"], log = T))
  }
  
  # return
  return(loglike)
}

r_logprior <- function(params, misc){
  ret <- 
    sum(dlnorm(params["group_sd"], 0, 5, log = TRUE)) +
    sum(dnorm(params[grepl("ccc", names(params))], 0, 10, log = TRUE)) +
    sum(dgamma2(params["dur_die"], 2, 1, log = TRUE)) +
    sum(dgamma2(params["dur_recover"], 30, 500, log = TRUE)) +
    sum(dgamma2(params[grepl("hosp", names(params))], 1, 0.7, log = TRUE))
  return(ret)
}
