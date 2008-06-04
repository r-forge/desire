##
## minimum_di.R - minimum desirability index
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <olafm@statistik.uni-dortmund.de>
##

minimumDI <- function(f, ...)
  UseMethod("minimumDI", f)

## Vector input
minimumDI.numeric <- function(f, ...) 
  min(f)

## Matrix input
minimumDI.matrix <- function(f, margin=1, ...)
  apply(f, margin, min)

## Array input
minimumDI.array <- function(f, margin=1, ...)
  apply(f, margin, min)

  
minimumDI.desire.function <- function(f, ...) {
  ev <- function(x)
    min(sapply(i, function(k) dfs[[k]](x[k])))

  dfs <- list(f, ...)
  if (!all(sapply(dfs, is.desirability)))
    stop("Not all supplied arguments are desirability functions.")
  
  i <- 1:length(dfs)
  class(ev) <- "desire.index"
  return(ev)
}

minimumDI.composite.desire.function <- function(f, ...) {
  ev <- function(x)
    min(sapply(dfs, function(f) f(x)))

  dfs <- list(f, ...)
  if (!all(sapply(dfs, is.composite.desirability)))
    stop("Not all supplied arguments are composite desirability functions.")
  
  class(ev) <- "desire.index"
  return(ev)
}
