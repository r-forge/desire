##
## stochastic.R - Stochastic (normaly distributed) desirability class
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <mersmann@statistik.uni-dortmund.de>
##

## Declare base clase 'desirability'
setClass("stochastic.desirability", representation("desirability",
                                                   d="desirability",
                                                   type="character",
                                                   sigma="numeric"))

setMethod("initialize", "stochastic.desirability", 
          function(.Object, ...) {
            .Object <- callNextMethod()
            if (.Object@sigma < 0)
              stop("'sigma' must be >= 0 but is ", sigma)
            return(.Object)
          })

setMethod("edesirability", "stochastic.desirability",
          function(dfun, y) {
            z <- sapply(y, function(yy) {
              mean(rdesirability(10000, dfun@d, yy, dfun@sigma))
            })
            return(z)
          })

sdesire <- function(d, sigma) {  
  val <- new("stochastic.desirability",
             type="Stochastic desirability",
             d=d, sigma=sigma)
  return(val)
}
