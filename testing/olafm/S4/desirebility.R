##
## desirability.R - Base class of all desirability functions
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <mersmann@statistik.uni-dortmund.de>
##

## Declare base clase 'desirability'
setClass("desirability", representation("function", type="character"))

##
## desirability.initialize - Constructor
##
## Implementation:
##
##  'desirability' is a subclass of 'function', as such its .Data slot should contain
##  a function definition. Here we initialize it to a wrapper which calls the corresponding
##  edesirability method. This is mainly a convenience for the user, so (s)he can type
##    d <- some.desirability()
##    ## evaluate d in 2
##    d(2) 
##    ## alternatevly
##    edesirability(d, 2)
##
setMethod("initialize", "desirability", 
          function(.Object, ...) {
            .Object <- callNextMethod()
            .Object@.Data <- function(value) { return(edesirability(value, sys.function(-1))) }
            return(.Object)
          })

## Return interesting range of desirability
##
## There is no sane way to find an interesting range. Just return a default one.
setMethod("irange", "desirability",
          function(f) {
            return(c(0, 10))
          })

setMethod("plot", "desirability",
          function(x, y, ..., xlim=irange(x), ylim=c(0, 1), xlab="Value", ylab="desirability", n=501) {
            plot.new()
            plot.window(xlim, ylim)
            box() ; axis(1) ; axis(2)
            title(main=paste("Type: ", x@type), xlab=xlab, ylab=ylab)
            xrng <- par("usr")[1:2]
            p <- seq(xrng[1], xrng[2], length.out=n)
            y <- x(p)
            abline(h=c(0, 1), col="grey", lty=2)
            lines(p, y, ...)
          })

