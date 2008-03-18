##
## harrington1.R - One sided Harrington type desirabilities
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <mersmann@statistik.uni-dortmund.de>
##

## One sided Harrington desirability
setClass("desirability.harrington.1",
         representation("desirability",
                        y1="numeric", d1="numeric",
                        y2="numeric", d2="numeric",
                        b0="numeric", b1="numeric"))

setMethod("initialize", "desirability.harrington.1",
          function(.Object, ...) {
            .Object <- callNextMethod()
            ## Solve for constants b0, b1
            ## See Heike Trautmanns Diss. p. 13-14
            X <- cbind(1, c(.Object@y1, .Object@y2))
            b <- solve(X, -log(-log(c(.Object@d1, .Object@d2))))
            .Object@b0 <- b[1]; .Object@ b1 <- b[2]
            .Object@type <- "One sided Harrington"
            return(.Object)
          })

## 
## irange - return range so desirability will be plotted from 0.05 to 0.95.
##
setMethod("irange", "desirability.harrington.1",
          function(f) {
            l <-  uniroot(function(x) f(x) - 0.05,
                          c(-1e10, 1e10),
                          tol=.Machine$double.eps)$root
            u <-  uniroot(function(x) f(x) - 0.95,
                          c(-1e10, 1e10),
                          tol=.Machine$double.eps)$root
            return(sort(c(l, u)))
          })

setMethod("show", "desirability.harrington.1",
          function(object) {
            y <- format(c(object@y1, object@y2))
            d <- format(c(object@d1, object@d2))
            cat("One sided Harrington desirability", "\n",
                " Parameters:", "\n",
                "  y1 = ", y[1], "   d1 = ", d[1], "\n",
                "  y2 = ", y[2], "   d2 = ", d[2], "\n", sep="")
            invisible(NULL)
          })

setMethod("edesirability", "desirability.harrington.1",
          function(y, dfun) {
            ys <- dfun@b0 + dfun@b1*y
            return(exp(-exp(-ys)))
          })

setMethod("ddesirability", c("numeric", "desirability.harrington.1", "numeric", "numeric"),
          function(x, dfun, mean, sd) {
              mu.t <- -(dfun@b0 + dfun@b1*mean)
              sq.t <- (dfun@b1 * sd)^2              
              ## c.f. Trautmann Diss p. 51
              c <- -(sqrt(2 * pi * sq.t) * log(x) * x)^(-1)
              ee <- -(2*sq.t)^(-1) * ((log(-log(x)) - mu.t)^2)
              return(c * exp(ee))
          })

setMethod("pdesirability", c("numeric", "desirability.harrington.1", "numeric", "numeric"),
          function(q, dfun, mean, sd) {
            mu.t <- -(dfun@b0 + dfun@b1*mean)
            sd.t <- dfun@b1 * sd
            ## c.f. Trautmann Diss p. 54
            return(1 - pnorm((log(-log(q)) - mu.t)/sd.t))
          })

setMethod("qdesirability", c("numeric", "desirability.harrington.1", "numeric", "numeric"),
          function(p, dfun, mean, sd) {
            mu.t <- -(dfun@b0 + dfun@b1*mean)
            sd.t <- dfun@b1 * sd
            z <- qnorm(1-p)            
            ## c.f. Trautmann Diss p. 54
            return(exp(-exp(sd.t * z + mu.t)))
          })

harrington1 <- function(y1, d1, y2, d2) {  
  val <- new("desirability.harrington.1",
             y1=y1, d1=d1, y2=y2, d2=d2)
  return(val)
}
