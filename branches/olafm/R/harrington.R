##
## harrington.R - Harrington type desiraility functions
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##

harrington1 <- function(y1, d1, y2, d2) {
  ev <- function(y) {
    ys <- b0 + b1*y
    return(exp(-exp(-ys)))
  }
  ## Solve for constants b0, b1
  ## See Heike Trautmanns Diss. p. 13-14
  X <- cbind(1, c(y1, y2))
  b <- solve(X, -log(-log(c(d1, d2))))
  b0 <- b[1];  b1 <- b[2]
  
  class(ev) <- c("desire.function.harrington1", "desire.function")
  attr(ev, "desire.type") <- "Two sided Harrington"
  attr(ev, "y.range") <- c(y1, y2)
  ## Remove cruft to save space
  rm(b, X)
  return(ev)
}

harrington2 <- function(LSL, USL, n) {
  ev <- function(y) {
    ys <- (2*y - s)/d
    return(exp(-abs(ys)^n))
  }
  if (USL <= LSL)
    stop("LSL must be strictly less than USL.")
         
  ## Precompute:
  s <- USL + LSL
  d <- USL - LSL
  class(ev) <- c("desire.function.harrington2", "desire.function")
  attr(ev, "desire.type") <- "Two sided Harrington"
  attr(ev, "y.range") <- c(LSL-d, USL + d)
  return(ev)
}

plot.desire.function.harrington2 <- function(f, ...) {
  e <- environment(f)
  lsl <- get("LSL", envir=e)
  usl <- get("USL", envir=e)
  plot.desire.function(f, ...)
  abline(v=c(lsl, usl), col="grey", lty=2)
  mtext(c("LSL", "USL"), at=c(lsl, usl), line=.5)
}


ddesire.desire.function.harrington1 <- function(x, f, mean=0, sd=1) {
  e <- environment(f)
  b0 <- get("b0", envir=e)
  b1 <- get("b1", envir=e)

  mu.t <- -(b0 + b1*mean)
  sq.t <- b1^2*sd^2
    
  ## c.f. Trautmann Diss p. 51
  c <- -(sqrt(2*pi*sq.t) * log(x) * x)^(-1)
  ee <- -(2*sq.t)^(-1) * ((log(-log(x)) - mu.t)^2)
  return(c * exp(ee))
}

pdesire.desire.function.harrington1 <- function(q, f, mean=0, sd=1) {
  e <- environment(f)
  b0 <- get("b0", envir=e)
  b1 <- get("b1", envir=e)
  mu.t <- -(b0 + b1*mean)
  sd.t <- b1*sd
  
  ## c.f. Trautmann Diss p. 54
  return(1 - pnorm((log(-log(q)) - mu.t)/sd.t))
}

qdesire.desire.function.harrington1 <- function(p, f, mean=0, sd=1) {
  e <- environment(f)
  b0 <- get("b0", envir=e)
  b1 <- get("b1", envir=e)
  mu.t <- -(b0 + b1*mean)
  sd.t <- b1*sd
  z <- qnorm(1-p)
  
  ## c.f. Trautmann Diss p. 54
  return(exp(-exp(sd.t * z + mu.t)))
}

h <- harrington1(1, .2, 6, .5)
plot(density(r1(1000, h)))
lines(density(r2(1000, h)), col="red")

rdesire.desire.function.harington1 <- function(

