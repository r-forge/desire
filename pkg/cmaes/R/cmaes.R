##
## cma.es - covariance matrix adapting evolutionary strategy
##
## This function implements an evolutionary strategy with covariance
## matrix adaption for function optimiziation. It is based on Matlab
## code by Nikolas Hansen. The call interface is modeled after the
## optim() function in the stats package.
##
## Input:
##   par : Inital solution
##   fn : function to minimize
##   ... : Additional arguments passed to fun()
##   control : list with options
##
##   Possible members of the control list are:
##
##   fnscale : Minimize fun() * fnscale
##   stopfitness : Minimum decrease in fitness before algorithm terminates
##   maxit : Maximum number of iterations / generations
##   sigma : Inital variance estimates
##
## Output:
##   A list with the following members:
##
##   par : Best parameter settings found
##   value : Fitness function value for par
##   counts : Number of function evaluations
##   convergence: 1 if the algorithm stopped due to meeting the
##     desiresd fitness value, 0 else.
##
##   For a detailed description see ?optim.
##
## Author:
##  Olaf Mersmann (OME) <olafm@statistik.tu-dortmund.de>
##
## Changes:
##  2009-03-18: (OME)
##   * Ported Matlab version by N. Hansen to R
##
##  2009-04-17: (OME)
##   * Code cleanup
##   * Better comments / documentation
##
## 2009-04-24: (OME)
##   * Turn file into R package

cma.es <- function(par, fn, ..., lower, upper, control=list()) {
  norm <- function(x)
    sqrt(crossprod(x))

  controlParam <- function(name, default) {
    v <- control[[name]]
    if (is.null(v))
      return (default)
    else
      return (v)
  }

  ## Inital solution:
  xmean <- par
  N <- length(xmean)

  ## Box constraints:
  if (missing(lower))
    lower <- rep(-Inf, N)
  else if (length(lower) == 1)  
    lower <- rep(lower, N)

  if (missing(upper))
    upper <- rep(Inf, N)
  else if (length(upper) == 1)  
    upper <- rep(upper, N)

  ## Parameters:
  trace       <- controlParam("trace", FALSE)  
  fnscale     <- controlParam("fnscale", 1)
  stopfitness <- controlParam("stopfitness", -Inf)
  maxiter     <- controlParam("maxit", 100 * N^2)
  sigma       <- controlParam("sigma", rep(0.5, N))
  if (length(sigma) == 1)
    sigma <- rep(sigma, N)

  ## Strategy parameter setting: Selection  
  lambda  <- controlParam("lambda", 4+floor(3*log(N)))
  mu      <- controlParam("mu", floor(lambda/2))
  weights <- controlParam("weights", log(mu+1) - log(1:mu))
  weights <- weights/sum(weights)
  mueff   <- controlParam("mueff", sum(weights)^2/sum(weights^2))
  cc      <- controlParam("ccum", 4/(N+4))
  cs      <- controlParam("cs", (mueff+2)/(N+mueff+3))
  mucov   <- controlParam("ccov.mu", mueff)
  ccov    <- controlParam("ccov.1",
                        (1/mucov) * 2/(N+1.4)^2
                        + (1-1/mucov) * ((2*mucov-1)/((N+2)^2+2*mucov)))
  damps <- controlParam("damps",
                        1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs)

  ## Safty checks:
  stopifnot(length(upper) == N)  
  stopifnot(length(lower) == N)
  stopifnot(all(lower <= upper))
  stopifnot(length(sigma) == N)
  
  ## Initialize dynamic (internal) strategy parameters and constants
  pc <- rep(0.0, N)
  ps <- rep(0.0, N)
  B <- diag(N)
  D <- diag(N)
  BD <- B %*% D
  C <- BD %*% t(BD)

  chiN <- sqrt(N) * (1-1/(4*N)+1/(21*N^2))
  
  iter <- 0L      ## Number of iterations
  counteval <- 0L ## Number of function evaluations

  ## Preallocate work arrays:
  arx <- matrix(0.0, nrow=N, ncol=lambda)
  arfitness <- numeric(lambda)
  while (iter < maxiter) {
    iter <- iter + 1L
    
    arz <- matrix(rnorm(N*lambda), nrow=N, ncol=lambda)
    for (k in 1:lambda) {
      ## Transform uncorrelated N(0,1) random vectors to the desired
      ## N(xmean, sigma) distribution:
      arx[,k] <- xmean + sigma * (BD %*% arz[,k])
      ## Enforce bounds:
      lidx <- arx[,k] < lower
      if (any(lidx))
        arx[lidx, k] <- lower[lidx]
      uidx <- arx[,k] > upper
      if (any(uidx))
        arx[uidx, k] <- upper[uidx]
      ## Calculate fitness:
      arfitness[k] <- fn(arx[,k], ...) / fnscale
      counteval <- counteval + 1L;
    }
    
    ## Order fitness:
    arindex <- order(arfitness)
    arfitness <- arfitness[arindex]

    aripop <- arindex[1:mu]
    selx <- arx[,aripop]
    xmean <- drop(selx %*% weights)
    selz <- arz[,aripop]
    zmean <- drop(selz %*% weights)

    ## Cumulation: Update evolutionary paths
    ps <- (1-cs)*ps + sqrt(cs*(2-cs)*mueff) * (B %*% zmean)
    hsig <- drop((norm(ps)/sqrt(1-(1-cs)^(2*counteval/lambda))/chiN) < (1.4 + 2/(N+1)))
    pc <- (1-cc)*pc + hsig * sqrt(cc*(2-cc)*mueff) * drop(BD %*% zmean)

    ## Adapt Covariance Matrix:
    BDz <- BD %*% selz
    C <- (1-ccov) * C + ccov * (1/mucov) *
      (pc %o% pc + (1-hsig) * cc*(2-cc) * C) +
        ccov * (1-1/mucov) * BDz %*% diag(weights) %*% t(BDz)
    
    ## Adapt step size sigma:
    sigma <- sigma * exp((norm(ps)/chiN - 1)*cs/damps)
    
    e <- eigen(C, symmetric=TRUE)
    B <- e$vectors
    D <- if (length(e$values) > 1L)
      diag(sqrt(e$values))
    else
      as.matrix(sqrt(e$values))
    BD <- B %*% D

    ## break if fit:
    if (arfitness[1] <= stopfitness)
      break;

    ## Escape from flat-land:
    if (arfitness[1] == arfitness[min(1+floor(lambda/2), 2+ceiling(lambda/4))]) { 
      sigma <- sigma * exp(0.2+cs/damps); 
      warning("Flat fitness function. Increasing sigma.")
    }
  }
  cnt <- vector("integer", 2L)
  cnt[1] <- as.integer(counteval)
  cnt[2] <- NA
  names(cnt) <- c("function", "gradient")

  ## Currently only the 'best' solution is returned, but since the
  ## CMA-ES is not a 'global' algorithm in the sense SANN is, all
  ## other solutions are supposed to be 'near' to the best on in the
  ## parameter space.
  res <- list(par=arx[, arindex[1]],
              value=arfitness[1] * fnscale,
              counts=cnt,
              convergence=ifelse(iter >= maxiter, 1L, 0L),
              message=NULL)
  return(res)
}
