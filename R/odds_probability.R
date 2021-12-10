# Probability -> odds
p <- function(odds){
  odds / (1 + odds)
}

# Odds -> probability
o <- function(prob){
  prob / (1 - prob)
}