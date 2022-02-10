# Gamma re-parameterised with mean and variance.

rgamma2 <- function(n, mean, var){
  shape = (mean ^ 2) / var
  scale = var / mean
  rgamma(n, shape = shape, scale = scale)
}

dgamma2 <- function(x, mean, var, log = FALSE){
  shape = (mean ^ 2) / var
  scale = var / mean
  dgamma(x, shape = shape, scale = scale, log = log)
}