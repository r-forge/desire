##
## desirability.R - R/P/Q/D functions for desirability functions
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##

## FIXME: N(0,1) realistic? Do we want defaults at all?
pdesire <- function(q, f, mean=0, sd=1) {
  UseMethod("pdesire", f)
}

qdesire <- function(p, f, mean=0, sd=1) {
  UseMethod("qdesire", f)
}

ddesire <- function(x, f, mean=0, sd=1) {
  UseMethod("ddesire", f)
}

rdesire <- function(n, f, mean=0, sd=1) {
  UseMethod("ddesire", f)
}

## Is this 'sane' or should we use the inversion rule?
rdesire.default <- function(n, f, mean=0, sd=1) {
  return(f(rnorm(n, mean, sd)))
}

## Alternative:
## rdesire.default <- function(n, f, mean=0, sd=1) {
##   return(qdesire(runif(n), f, mean, sd))
## }
