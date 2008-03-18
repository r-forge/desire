##
## harrington2.R - Two sided Harrington type desirabilites
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <mersmann@statistik.uni-dortmund.de>
##


## Two sided Harrington desirability
setClass("desirability.harrington.2",
         representation("desirability",
                        LSL="numeric", USL="numeric", n="numeric"))

setMethod("initialize", "desirability.harrington.2",
          function(.Object, ...) {
            .Object <- callNextMethod()
            ## Check for sane USL / LSL
            if (.Object@USL <= .Object@LSL)
              stop("LSL must be strictly less than USL.")
            return(.Object)
          })

setMethod("irange", "desirability.harrington.2",
          function(f) {
            s <- f@USL + f@LSL
            d <- f@USL - f@LSL
            z <- d*(-log(0.05))^(1/f@n)
            l <- .5 * (s - z)
            u <- .5 * (s + z)
            return(c(l, u))
          })

setMethod("show", "desirability.harrington.2",
          function(object) {
            s <- format(c(object@LSL, object@USL))
            cat("Two sided Harrington desirability", "\n",
                " Parameters:", "\n",
                "  LSL = ", s[1], "   USL = ", s[2], "\n", sep="")
            invisible(NULL)
          })

setMethod("edesirability", "desirability.harrington.2",
          function(y, dfun) {
            s <- dfun@USL + dfun@LSL
            d <- dfun@USL - dfun@LSL
            ys <- (2*y - s)/d
            return(exp(-abs(ys)^dfun@n))
          })

setMethod("ddesirability", c("numeric", "desirability.harrington.2", "numeric", "numeric"),
          function(x, dfun, mean, sd) {
            ## See Trautmann p. 51
            s <- dfun@USL + dfun@LSL
            d <- dfun@USL - dfun@LSL
            n <- dfun@n ; n0 <- 1/n ; n1 <- n0 - 1

            mu.t <- (2/d) * mean - s / d
            sd.t <- (2/d) * sd

            c0 <- sqrt(2*pi) * sd.t * x * n
            c1 <- -log(x); c1a <- c1^n0; c1b <- c1^n1            
            c2 <- 2*sd.t^2
            c3 <- (c1a - mu.t)^2 / c2
            c4 <- (c1a + mu.t)^2 / c2

            r <- c1b/c0 * (exp(-c3) + exp(-c4))
            return(ifelse(x < 0, 0, ifelse(x > 1, 0, r)))
          })

setMethod("pdesirability", c("numeric", "desirability.harrington.2", "numeric", "numeric"),
          function(q, dfun, mean, sd) {
            ## See Trautmann p. 51
            s <- dfun@USL + dfun@LSL
            d <- dfun@USL - dfun@LSL
            n <- dfun@n
            mu.t <- 2 / d * mean - s / d
            sd.t <- (2 / d) * sd

            ## FIXME: 1/n - 1 or 1/(n-1)
            n0 <- 1/n
            c1 <- -log(x)^n0
            
            return(2 - pnorm((c1 - mu.t)/sd.t) - pnorm((c1 + mu.t)/sd.t))
          })

harrington2 <- function(LSL, USL, n) {
  val <- new("desirability.harrington.2",
             type="Two sided Harrington",
             LSL=LSL, USL=USL, n=n)
  return(val)
}

