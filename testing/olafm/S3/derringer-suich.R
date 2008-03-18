##
## derringer-suich.R - Derringer-Suich type desirability functions
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##

derringerSuich <- function(y, d, beta) {
  ev <- function(t) {
    i <- findInterval(t, yx, rightmost.closed=TRUE)
    y.i <- yx[i+1] ; y.im <- yx[i]
    d.i <- dx[i+1] ; d.im <- dx[i]
    b.i <- bx[i]
    ifelse(d.im <= d.i,
           d.im + (d.i  - d.im)*((t - y.im)/(y.i  - y.im))^b.i,
           d.i  + (d.im - d.i )*((t - y.i )/(y.im - y.i ))^b.i)
  }
  
  n <- length(y)
  if (length(d) != n)
    stop("Number of desirabilities does not match number of data points.")
  if (length(beta) != (n-1))
    stop("Number of weights does not match number of data points.")
  if (is.unsorted(y))
    stop("Data points 'y' not ordered.")

  if (any(d < 0 | d > 1))
    stop("Not all desirabilities in the range 0 to 1,")
  if (any(beta <= 0))
    stop("Not all weights are positive.")
  
  ## Extend range to eliminate < y_min and > y_max conditions:
  yx <- c(-Inf, y, +Inf)
  dx <- c(0, d, 0)
  bx <- c(0, beta, 1)

  class(ev) <- "desire.function"
  attr(ev, "desire.type") <- "Derringer-Suich"
  attr(ev, "y.range") <- range(y[is.finite(y)])
  ## Remove unnecessary variables, since they will be saved in ev's environment. 
  rm(y, d, beta, n)
  return(ev)
}

