# Probability from odds
p <- function(odds){
  odds / (1 + odds)
}

# Odds from probability
o <- function(prob){
  prob / (1 - prob)
}