##
## loglognorm.R - Interface to loglognorm.c
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <olafm@statistik.uni-dortmund.de>
##

dloglognorm <- function(x, mean=0, sd=1)
  .External("dloglognorm", x, mean, sd, PACKAGE="loglognorm")

ploglognorm <- function(q, mean=0, sd=1)
  .External("ploglognorm", q, mean, sd, PACKAGE="loglognorm")

qloglognorm <- function(p, mean=0, sd=1)
  .External("qloglognorm", p, mean, sd, PACKAGE="loglognorm")

## FIXME: Faster way?
rloglognorm <- function(n, mean=0, sd=1)
  .External("qloglognorm", runif(n), mean=mean, sd=sd, PACKAGE="loglognorm")

mloglognorm <- function(mean, sd, moment)  
  .External("mloglognorm", mean, sd, moment, PACKAGE="loglognorm")

eloglognorm <- function(mean, sd)
  .External("mloglognorm", mean, sd, rep(1, length(mean)), PACKAGE="loglognorm")

vloglognorm <- function(mean, sd) {
  m1 <- mloglognorm(mean, sd, rep(1, length(mean)))
  m2 <- mloglognorm(mean, sd, rep(2, length(mean)))
  return (m2 - m1^2)
}

