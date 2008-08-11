set.seed(42)

require(mco)

normalizeFront <- function(front, min, max)
  t((t(front) - min)/(max - min))

## Squared distance
distance2 <- function(x, y) {
  (x - y) %*% (x - y)
}

## Squared distance to front
distanceToFront2 <- function(x, front) {
  min(sapply(1:nrow(front), function(i) { z <- (x - front[i,]); z %*% z }))
}

generationalDistance <- function(x, o) {
  front <- x$value[x$pareto.optimal,]
  if (is.matrix(o)) {
    o <- truefront
  } else if ("nsga2" %in% class(o)) {
    truefront <- o$value[o$pareto.optimal,]
  } else {
    stop("Don't know how to interpret true front of class ", class(o), ".")
  }

  ## Normalize front:
  maxval <- apply(truefront, 2, max)
  minval <- apply(truefront, 2, min)

  nfront <- normalizeFront(front, minval, maxval)
  ntruefront <- normalizeFront(truefront, minval, maxval)

  ## Calculate criterion:
  d <- sapply(1:nrow(nfront), function(i) distanceToFront2(nfront[i,], ntruefront))
  return(sqrt(sum(d))/nrow(nfront))
}

generalizedSpread <- function(x, o) {
  front <- x$value[x$pareto.optimal,]
  if (is.matrix(o)) {
    o <- truefront
  } else if ("nsga2" %in% class(o)) {
    truefront <- o$value[o$pareto.optimal,]
  } else {
    stop("Don't know how to interpret true front of class ", class(o), ".")
  }
  
  ## Normalize front:
  maxval <- apply(truefront, 2, max)
  minval <- apply(truefront, 2, min)

  nfront <- normalizeFront(front, minval, maxval)
  ntruefront <- normalizeFront(truefront, minval, maxval)

  K <- nrow(nfront)
  N <- nrow(ntruefront)
  ## Calculate extreme values:
  nobj <- ncol(front)
  
  extreme <- matrix(0, ncol=nobj, nrow=nobj)
  for (i in 1:nobj) {
    o <- order(ntruefront[,i])
    for (j in 1:nobj) {
      extreme[i,j] <- ntruefront[o,][N, j]
    }
  }
  ## Sort 'true' pareto front:
  for (i in nobj:1)
    nfront <- nfront[order(nfront[,i]),]

  if (distance2(nfront[1,], nfront[N,]) == 0) {
    return (0.0)
  } else {
    dmean <- mean(sapply(1:K, function(i) sqrt(distanceToFront2(nfront[i,], nfront[-i,]))))
    dextr <- sum(sapply(1:nobj, function(i) sqrt(distanceToFront2(extreme[i,], nfront))))
    mean <-  sum(sapply(1:K, function(i) sqrt(distanceToFront2(nfront[i,], nfront[-i,]))-dmean))
    return ((dextr + mean)/(dextr + K*dmean))                 
  }
}

################################################################################
## Example:
sch1 <- function(x, ...) {
  c(x[1]^2, (x[1] - 2)^2)
}

res1 <- nsga2(sch1, 1, 2,
             lower.bounds=0, upper.bounds=1,
              popsize=20, generations=10)

res2 <- nsga2(sch1, 2, 2,
              lower.bounds=0, upper.bounds=1,
              popsize=16, generations=seq(2, 100, by=5))
## plot(res2)


n <- length(res2)
GD <- sapply(1:n, function(x) generationalDistance(res2[[x]], res2[[n]]))
GS <- sapply(1:n, function(x) generalizedSpread(res2[[x]], res2[[n]]))


 
