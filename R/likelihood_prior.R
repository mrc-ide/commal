r_loglike <- function(params, data, misc) {
  country_capacity <- params[grepl("ccc", names(params))]
  country_hosp <-params[grepl("hosp_", names(params))]
  block <- unlist(misc["block"])
  
  if(block %in% 1:22){
    prob_dhs <- misc$model_function(pfpr = data$dhs[[block]]$pfpr,
                                    distance = data$dhs[[block]]$distance,
                                    global_capacity = params["global_capacity"],
                                    country_capacity = country_capacity[block],
                                    distance_beta = params["distance_beta"],
                                    pfpr_beta = params["pfpr_beta"],
                                    shift = params["shift"])
    loglike <- sum(dbinom(data$dhs[[block]]$sma, 1, prob_dhs, log = T))
  }
  
  if(block %in% 23:25){
    paton_block <- block - 22
    prob_paton <- misc$model_function(pfpr = data$paton[[paton_block]]$pfpr,
                                      distance = data$paton[[paton_block]]$distance,
                                      global_capacity = params["global_capacity"],
                                      country_capacity = country_capacity[paton_block],
                                      distance_beta = params["distance_beta"],
                                      pfpr_beta = params["pfpr_beta"],
                                      shift = params["shift"])
    prob_paton <- misc$sma_prev_age_standardise(prob_paton)
    durations <- params[c("dur_recover", "dur_tx", "dur_die")]
    weights <- cbind(
      p_recover = (1 - data$paton[[paton_block]]$act) * (1 - params["cfr"]),
      p_tx = data$paton[[paton_block]]$act,
      p_die = (1 - data$paton[[paton_block]]$act) * params["cfr"]
    )
    mean_duration <- apply(weights, 1, weighted.mean, x = durations)
    inc <- data$paton[[paton_block]]$py * 365 * inc1(prob_paton, 1 / mean_duration)
    hosp_inc <- inc / country_hosp[paton_block]
    loglike <- sum(dpois(data$paton[[paton_block]]$sma, hosp_inc, log = T))
  }
  
  if(block == 26){
    loglike <- sum(dnorm(country_capacity, 0, params["group_sd"], log = T))
  }
  
  # return
  return(loglike)
}

r_logprior <- function(params, misc){
  ret <- 
    sum(dlnorm(params["group_sd"], 0, 5, log = TRUE)) +
    sum(dnorm(params[grepl("ccc", names(params))], 0, 10, log = TRUE)) +
    sum(dgamma2(params["dur_recover"], 315, 1000, log = TRUE)) +
    sum(dgamma2(params["dur_tx"], 20, 100, log = TRUE)) +
    sum(dgamma2(params["dur_die"], 20, 100, log = TRUE)) +
    sum(dbeta(params["cfr"], 5, 5, log = TRUE))# +
  #sum(dlnorm(params[grepl("hosp", names(params))], 0, 1, log = TRUE))
  return(ret)
}
