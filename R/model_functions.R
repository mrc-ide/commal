# Estimate prob of SMA
gompertz <- function(pfpr, global_capacity, country_capacity, pfpr_beta, shift){
  carrying_capacity <- rlogit(global_capacity + country_capacity)
  est <- carrying_capacity * exp(-exp(shift - pfpr_beta * pfpr))
  return(est)
}

# Estimate the probability from log odds
rlogit <- function(log_odds){
  1 / (1 + exp(-log_odds))
}

dic <- function(loglikelihood){
  d <- -2 * loglikelihood
  0.5 * var(d) + mean(d)
}