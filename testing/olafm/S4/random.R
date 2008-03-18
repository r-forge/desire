##
## random.R - Default generics for stochastic desirabilities
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <mersmann@statistik.uni-dortmund.de>
##

setMethod("ddesirability",
          c("numeric", "desirability", "numeric", "numeric"),
          function(x, dfun, mean, sd) {
            n <- 200000
            warning("Using ", n, " sample approximation for desirability quantile")
            z <- rdesirability(n, dfun, mean, sd)
            ## FIXME: Ugly hack, because density() does not let me choose evaluation points
            ##
            ## Use kernel estimate of density
            d <- density(z, from=0, to=1, n=2048)
            ## And linearly interpolate 
            f <- approxfun(d$x, d$y)
            return(f(x))
          })

##
## pdesirability - CDF of normal desirability
##
## Implementation:
##
##   The CDF is approximated by the empirical CDF of 10000 random samples from the
##   desirabilities distribution. This is neither stable nor elegant, but the only 
##   doable thing, given that we cannot assume to have a density which we could 
##   numerically integrate.
##
setMethod("pdesirability", c("numeric", "desirability", "numeric", "numeric"),
          function(q, dfun, mean, sd) {
            warning("Using 10000 sample approximation for desirability quantile")
            x <- rdesirability(10000, dfun, mean, sd)
            p <- ecdf(x)(q)
            return(p)
          })

## 
## qdesirability - Quantile function of normal desirability
##
## Implementation:
##   
##   The quantile function is approximated using the inverse of the empirical CDF
##   of 10000 random samples from the underlying distribution.
##
setMethod("qdesirability", c("numeric", "desirability", "numeric", "numeric"),
          function(p, dfun, mean, sd) {
            warning("Using 10000 sample approximation for desirability quantile")
            x <- rdesirability(10000, dfun, mean, sd)
            return(quantile(x, p))
          })

##
## rdesirability - Random numbers from a normal desirability
##
## Implementation:
##
##   'n' N(mean, sd) distributed random samples are pulled and transformed using 
##   the desirability function.
##
setMethod("rdesirability",
          c("numeric", "desirability", "numeric", "numeric"),
          function(n, dfun, mean, sd) {
            z <- rnorm(n, mean, sd)
            return(dfun(z))
          })
