##
## index.R -  desirebility index functions
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##

desirebilityIndex <- function(df, type=c("geometric", "maximin")) {
  ev <- function(...) {
    x <- c(...)
    d <- sapply(df, function(f) f(x))
    i <- apply(d, 1, prod)
    return(i)
  }  
  return(ev)
}
