##
## mean_di.R - weighted mean desirability index
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <olafm@statistik.uni-dortmund.de>
##

meanDI <- function(f, ..., weights=1) 
  UseMethod("meanDI", f)

## Vector input
meanDI.numeric <- function(f, ..., weights) {
  weights <- weights/sum(weights)
  mean(f*weights)
}
  
## Matrix input
meanDI.matrix <- function(f, margin=1, ..., weights=1) {
  weights <- weights/sum(weights)
  apply(f, margin, function(x) mean(x*weights))
}

## Array input
meanDI.array <- function(f, margin=1, ..., weights=1)  {
  weights <- weights/sum(weights)
  apply(f, margin, function(x) mean(x * weights))
}

meanDI.desire.function <- function(f, ..., weights=1) {
  weights <- weights/sum(weights)
  ev <- function(x)
    mean(sapply(i, function(k) dfs[[k]](x[k])) * weights)
  
  dfs <- list(f, ...)
  if (!all(sapply(dfs, is.desirability)))
    stop("Not all supplied arguments are desirability functions.")
  
  i <- 1:length(dfs)
  class(ev) <- "desire.index"
  return(ev)
}

meanDI.composite.desire.function <- function(f, ..., weights=1) {
  weights <- weights/sum(weights)
  ev <- function(x)
    mean(sapply(dfs, function(f) f(x)) * weights)
  
  dfs <- list(f, ...)
  if (!all(sapply(dfs, is.composite.desirability)))
    stop("Not all supplied arguments are composite desirability functions.")
  
  class(ev) <- "desire.index"
  return(ev)
}
