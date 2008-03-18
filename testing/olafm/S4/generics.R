##
## generics.r - Declare new generic functions
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <mersmann@statistik.uni-dortmund.de>
##

## irange - interesting range of a function
setGeneric("irange", function(f) { standardGeneric("irange") })

## desirability:
setGeneric("edesirability", function(y, dfun) { standardGeneric("edesirability") })
setGeneric("ddesirability", function(x, dfun, mean, sd) { standardGeneric("ddesirability") })
setGeneric("pdesirability", function(q, dfun, mean, sd) { standardGeneric("pdesirability") })
setGeneric("qdesirability", function(p, dfun, mean, sd) { standardGeneric("qdesirability") })
setGeneric("rdesirability", function(n, dfun, mean, sd) { standardGeneric("rdesirability") })
