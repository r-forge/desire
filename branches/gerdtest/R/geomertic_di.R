##
## geometric_di.R - geometric desirability index
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <olafm@statistik.uni-dortmund.de>
##

geometricDI <- function(f, ..., weights)
  UseMethod("geometricDI", f)

## Vector
geometricDI.numeric <- function(f, ..., weights) {
  if (missing(weights))
    weights <- rep(1, length(f))
  prod(f^weights)^(1/sum(weights))
}

## Matrix
geometricDI.matrix <- function(f, margin, ..., weights) {
  if (missing(weights))
    weights <- rep(1, dim(f)[margin])
  q <- 1 / sum(weights)
  apply(f, margin, function(x) prod(x^weights)^q)  
}

## Array
geometricDI.array <- function(f, margin, ..., weights) {
  if (missing(weights))
    weights <- rep(1, dim(f)[margin])
  q <- 1 / sum(weights)
  apply(f, margin, function(x) prod(x^weights)^q)  
}

## Plain desirability
geometricDI.desire.function <- function(f, ..., weights) {
  ev <- function(x) 
    prod(sapply(i, function(k) dfs[[k]](x[k])^weights[k]))^q

  dfs <- list(f, ...)
  if (!all(sapply(dfs, is.desirability)))
    stop("Not all supplied arguments are desirability functions.")
  
  i <- 1:length(dfs)
  if (missing(weights))
    weights <- rep(1, length(dfs))
  q <- 1/sum(weights)
  class(ev) <- "desire.index"
  return(ev)
}

geometricDI.composite.desire.function <- function(f, ..., weights) {
 ev <- function(x) {
   m <- sapply(i, function(k) dfs[[k]](x)^weights[k])
   ## Deal with lm etc. which might return a vector instead of a single value.   
   if (is.matrix(m)) {
     r <- apply(m, 1, prod)^q
     ## Copy rownames if input has any.
     if (!is.null(rownames(x))) 
       names(r) <- rownames(x)
   } else {
     r <- prod(m)^q
   }
   return(r)
 }

  dfs <- list(f, ...)
  if (!all(sapply(dfs, is.composite.desirability)))
    stop("Not all supplied arguments are composite desirability functions.")
  
  i <- 1:length(dfs)
  if (missing(weights))
    weights <- rep(1, length(dfs))
  q <- 1/sum(weights)
  class(ev) <- "composite.desire.index"
  return(ev)
}

