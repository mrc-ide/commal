# Parameterised functional fits form Paton et al (2021)
paton_sma <- function(pfpr, alpha = -13.51, beta = 7.2, gamma = 7.94, py = 1000){
  exp(alpha + beta / (1 + exp(-gamma * pfpr)) + log(py))
}

paton_rd <- function(pfpr, alpha = -8.37, beta = 2.14, py = 1000){
  exp(alpha + beta * pfpr + log(py))
}

paton_cm <- function(alpha = -8.57, py = 1000){
  exp(alpha + log(py))
}