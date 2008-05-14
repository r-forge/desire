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

eloglognorm <- function(mean=0, sd=1)
  .External("eloglognorm", mean, sd, PACKAGE="loglognorm")

