##
## truncnorm.R - Interface to truncnorm.c
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <olafm@statistik.uni-dortmund.de>
##

dtruncnorm <- function(x, a, b, mean=0, sd=1)
  .External("dtruncnorm", x, a, b, mean, sd)

ptruncnorm <- function(q, a, b, mean=0, sd=1)
  .External("ptruncnorm", q, a, b, mean, sd)

qtruncnorm <- function(p, a, b, mean=0, sd=1) {  
  ## Taken out of package 'msm'
  ret <- numeric(length(p))
  ret[p == 1] <- Inf
  ret[p == 0] <- -Inf
  ret[p < 0 | p > 1] <- NaN
  ind <- (p > 0 & p < 1)
  if (any(ind)) {
    hind <- seq(along = p)[ind]
    h <- function(y) {
      (ptruncnorm(y, a, b, mean, sd) - p)[hind[i]]
    }
    ptmp <- numeric(length(p[ind]))
    for (i in 1:length(p[ind])) {
      interval <- c(-1, 1)
      while (h(interval[1]) * h(interval[2]) >= 0) interval <- interval + 
        c(-1, 1) * 0.5 * (interval[2] - interval[1])
      ptmp[i] <- uniroot(h, interval, tol = .Machine$double.eps)$root
    }
    ret[ind] <- ptmp
  }
  return (ret)
}

rtruncnorm <- function(n, a, b, mean=0, sd=1)
  .External("rtruncnorm", n, a, b, mean, sd)

etruncnorm <- function(a, b, mean, sd)
  .External("etruncnorm", a, b, mean, sd)

vtruncnorm <- function(a, b, mean, sd) {
  na <- length(a)
  if (na > 1) {
    if (length(mean) == 1)
      mean <- rep(mean, na)
    if (length(sd) == 1)
      sd <- rep(sd, na)
  }
  .External("etruncnorm", a, b, mean, sd)
}

