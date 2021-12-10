r_loglike <- function(params, data, misc) {
  ret <- 0
  country_level <- params[grepl("ccc", names(params))]
  block <- unlist(misc["block"])
  
  if(block < 23){
    # Convert to probability
    prob <- misc$model_function(pfpr = data[[block]]$pfpr,
                                distance = data[[block]]$distance,
                                a = params["a"],
                                b = params["b"],
                                c = params["c"],
                                d = country_level[block],
                                e = params["e"])
    # calculate log-probability of data
    ret <- ret + sum(dbinom(data[[block]]$sma, 1, prob, log = T))
  } else {
    ret <- ret + sum(dnorm(country_level, 0, params["group_sd"], log = T))
  }
  # return
  return(ret)
}

r_logprior <- function(params, misc){
  ret <- 
    sum(dlnorm(params["group_sd"], 0, 5, log = TRUE)) +
    sum(dnorm(params[grepl("ccc", names(params))], 0, 10, log = TRUE))
  return(0)
}

# r_loglike <- function(params, data, misc) {
#   ### DHS GLM ##################################################################
#   # Estimate log odds
#   log_odds <- misc$model_function(model_data = data$dhs_data, params = params)
#   # Convert to probability
#   prob <- rlogit(log_odds)
#   # calculate log-probability of data
#   ret_glm <- dbinom2(data$dhs_data$sma, prob)
#   ##############################################################################
#   
#   ### Paton comparison #########################################################
#   # Predict SMA prevalence for Paton PfPr
#   log_odds_predict <- misc$model_function(model_data = data$paton_data, params = params)
#   prob_predict <- rlogit(log_odds_predict)
#   # Age standardise (convert prev 6m-5y to prev 6m-9y)
#   prob_predict <- misc$sma_prev_age_standardise(prob_predict)
#   
#   durations <- params[c("dur_recover", "dur_tx", "dur_die")]
#   weights <- cbind(
#     p_recover = (1 - data$paton_data$act) * (1 - params["cfr"]),
#     p_tx =    data$paton_data$act,
#     p_die = (1 - data$paton_data$act) * params["cfr"]
#   )
#   mean_duration <- apply(weights, 1, weighted.mean, x = durations)
#   
#   inc <- data$paton_data$py * 365 * inc1(prob_predict, 1 / mean_duration)
#   # Convert to hospitalised incidence
#   ors_params <- params[grepl("ors", names(params))]
#   hosp_inc <- inc / ors_params[data$paton_data$countryn]
#   # calculate log-probability of Paton data
#   ret_inc <- sum(dpois(data$paton_data$sma, hosp_inc, log = T))
#   ##############################################################################
# 
#   ret <- ret_glm + ret_inc
#   
# 
#   # return
#   return(ret)
#}

# r_logprior <- function(params, misc){
#   ret <- 
#     sum(dnorm(params[c("intercept", "beta", "gamma", "beta2", "country_2", "country_3")], 0, 10, log = TRUE)) +
#     sum(dgamma2(params["dur_recover"], 315, 1000, log = TRUE)) +
#     sum(dgamma2(params["dur_tx"], 20, 100, log = TRUE)) +
#     sum(dgamma2(params["dur_die"], 20, 100, log = TRUE)) +
#     sum(dbeta(params["cfr"], 5, 5, log = TRUE)) +
#     sum(dlnorm(params[c("ors1", "ors2", "ors3")], 0, 1, log = TRUE))
# 
#   return(ret)
# }

