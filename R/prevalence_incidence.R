#' Estimate the point prevalence from incidence and rate of recovery.
#' 
#' Assumes system at equilibrium and a constant rate of recovery.
#'
#' @param incidence Incidence
#' @param recovery_rate Recovery rare
#'
#' @return Point prevalence.
inc_to_prev <- function(incidence, recovery_rate){
  incidence / (incidence + recovery_rate)
}

#' Estimate the incidence from point prevalence.
#' 
#' Assumes system at equilibrium and a constant rate of recovery.
#'
#' @param prevalence Point prevalence
#' @param recovery_rate Recovery rare
#' @param py Person years
#'
#' @return Incidence rate.
prev_to_inc <- function(prevalence, recovery_rate, py = 1000){
  py * 365 * (prevalence * recovery_rate / (1 - prevalence))
}
