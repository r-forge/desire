##
## functions.R - Test functions
##
## Author:
##  Olaf Mersmann (OME) <olafm@statistik.tu-dortmund.de>
##

genShiftedRosenbrock <- function(o, bias) {
  f <- function(x) {
    z <- x - o + 1
    d <- length(z)
    hz <- z[1:(d-1)]
    tz <- z[2:d]
    s <- sum(100 * (hz^2 - tz)^2 + (hz - 1)^2)
    return(s + bias)
  }
  return(f)    
}
