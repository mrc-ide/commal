#' Estimate the period prevalence.
#' 
#' Assumes system at equilibrium and a constant rate of recovery.
#'
#' @param incidence Incidence rate
#' @param recovery_rate Recovery rare
#' @param period_prevalence_duration Duration of prevalence survey period
#' @param target Target incidence. Only used if numerically solving to estimate
#' incidence for an observed prevalence.
#'
#' @return Period prevalence.
period_prev <- function(incidence, recovery_rate, period_prevalence_duration, target = NULL){
  pp <- (incidence / (incidence + recovery_rate)) + ((recovery_rate / (incidence + recovery_rate)) * (1 - exp(-incidence * period_prevalence_duration)))
  if(!is.null(target)){
    pp <- abs(pp - target)
  }
  return(pp)
}

#' Estimate the point prevalence fro incidence and rate of recovery.
#' 
#' Assumes system at equilibrium and a constant rate of recovery.
#'
#' @param incidence Incidence
#' @param recovery_rate Recovery rare
#'
#' @return Point prevalence.
prev <- function(incidence, recovery_rate){
  incidence / (incidence + recovery_rate)
}

#' Estimate the incidence from point prevalence.
#' 
#' Assumes system at equilibrium and a constant rate of recovery.
#'
#' @param prevalence Point prevalence
#' @param recovery_rate Recovery rare
#'
#' @return Incidence rate.
inc1 <- function(prevalence, recovery_rate){
  prevalence * recovery_rate / (1 - prevalence)
}

#' Estimate the incidence from period prevalence.
#' 
#' Assumes system at equilibrium and a constant rate of recovery. Converges to
#' `inc1()` if `period_prevalence_duration = 0`.
#'
#' @param prevalence Point prevalence
#' @param recovery_rate Recovery rare
#' @param period_prevalence_duration Duration of prevalence survey period
#'
#' @return Incidence rate.
inc2 <- function(period_prevalence, recovery_rate, period_prevalence_duration){
  inc <- rep(NA, length(period_prevalence))
  
  for(i in seq_along(period_prevalence)){
    inc[i] <- 
      optimise(period_prev,
               interval = c(0, 2),
               tol = .Machine$double.eps,
               recovery_rate = recovery_rate,
               period_prevalence_duration = period_prevalence_duration,
               target = period_prevalence[i])$minimum
  }
  return(inc)
}
