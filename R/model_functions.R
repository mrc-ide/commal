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

# Standardise estimate of SMA prevalence for age range
# Default parameters are for medium transmission non-seasonal setting. See
# Fitting_SMA_age_distribution.R for more options

# prevalence - SMA prevalence
# distribution - age distribution function
# p1 - distribution parameter 1
# p2 - distribution parameter 2
# age_in_lower - lower end of age band for input prevalence (months)
# age_in_upper - upper end of age band for input prevalence (months)
# age_out_lower - lower end of age band for output prevalence (months)
# age_out_upper - upper end of age band for output prevalence (months)
sma_prev_age_standardise <- function(prevalence,
                                     dist = plnorm, p1 = 2.83, p2 = 0.79, 
                                     age_in_lower = 6, age_in_upper = 5 * 12,
                                     age_out_lower = 3, age_out_upper = 9 * 12){
  
  # Denominator adjustment (assuming same number of children in each age group)
  age_band_multiplier <- (age_out_upper - age_out_lower) / (age_in_upper - age_in_lower)
  # Numerator adjustment (assuming given age distribution of cases)
  cases_multiplier <- (dist(age_out_upper, p1, p2) - dist(age_out_lower, p1, p2)) / 
    (dist(age_in_upper, p1, p2) - dist(age_in_lower, p1, p2))
  
  prevalence_adjustment <- cases_multiplier / age_band_multiplier
  prevalence <- pmin(1, prevalence * prevalence_adjustment)
  return(prevalence)
}

dic <- function(loglikelihood){
  d <- -2 * loglikelihood
  0.5 * var(d) + mean(d)
}