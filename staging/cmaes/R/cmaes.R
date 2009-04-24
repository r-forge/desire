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

cma.es <- function(par, fn, ..., control=list()) {
  norm <- function(x)
    sqrt(crossprod(x))
  
  xmean <- par
  N <- length(xmean)

  fnscale <- control$fnscale
  if (is.null(fnscale))
    fnscale <- 1
  
  stopfitness <- control$stopfitness
  if (is.null(stopfitness))
    stopfitness <- -Inf
  
  stopeval <- control$maxit
  if (is.null(stopeval))
    stopeval <- 100 * N^2
  
  sigma <- control$sigma
  if (is.null(sigma))
    sigma <- rep(0.5, N)
  if (length(sigma) == 1)
    sigma <- rep(sigma, N)
  stopifnot(length(sigma)==N)

  ## Strategy parameter setting: Selection  
  lambda <- 4+floor(3*log(N));  
  mu <- floor(lambda/2);        
  weights <- log(mu+1)-log(1:mu)
  weights <- weights/sum(weights) # Normalize weights
  mueff <- sum(weights)^2/sum(weights^2)
  
  cc = 4/(N+4)
  cs = (mueff+2)/(N+mueff+3)
  mucov = mueff
  ccov = (1/mucov) * 2/(N+1.4)^2 + (1-1/mucov) * ((2*mucov-1)/((N+2)^2+2*mucov))
  damps = 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs
  
  ## Initialize dynamic (internal) strategy parameters and constants
  pc <- rep(0.0, N)
  ps <- rep(0.0, N)
  B <- diag(N)
  D <- diag(N)
  BD <- B %*% D
  C <- BD %*% t(BD)

  chiN <- sqrt(N) * (1-1/(4*N)+1/(21*N^2))
  counteval <- 0
  while (counteval < stopeval) {
    arx <- matrix(0, nrow=N, ncol=lambda)
    ## arz <- matrix(0, nrow=N, ncol=lambda)
    arz <- matrix(rnorm(N*lambda), nrow=N, ncol=lambda)
    arfitness <- numeric(lambda)
    for (k in 1:lambda) {
      ## arz[,k] <- rnorm(N)
      
      ## Transform uncorrelated N(0,1) random vectors to the desired
      ## N(xmean, sigma) distribution:
      arx[,k] <- xmean + sigma * (BD %*% arz[,k])
      arfitness[k] <- fn(arx[,k], ...) / fnscale
      counteval <- counteval+1;
    }
    
    ## Order fitness:
    arindex <- order(arfitness)
    arfitness <- arfitness[arindex]
    xmean <- drop(arx[,arindex[1:mu]] %*% weights)
    zmean <- drop(arz[,arindex[1:mu]] %*% weights)

    ## Cumulation: Update evolutionary paths
    ps <- (1-cs)*ps + sqrt(cs*(2-cs)*mueff) * (B %*% zmean)
    hsig <- drop((norm(ps)/sqrt(1-(1-cs)^(2*counteval/lambda))/chiN) < (1.4 + 2/(N+1)))
    pc <- (1-cc)*pc + hsig * sqrt(cc*(2-cc)*mueff) * drop(BD %*% zmean)

    ## Adapt C:
    BDz <- BD %*% arz[,arindex[1:mu]]
    C <- (1-ccov) * C + ccov * (1/mucov) *
      ((pc %*% t(pc)) + (1-hsig) * cc*(2-cc) * C) +
        ccov * (1-1/mucov) * BDz %*% diag(weights) %*% t(BDz)
    
    ## Adapt step size sigma:
    sigma <- sigma * exp((cs/damps)*(norm(ps)/chiN - 1))
    
    e <- eigen(C, symmetric=TRUE)
    B <- e$vectors
    D <- if (length(e$values) > 1)
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
  cnt <- vector("integer", 2)
  cnt[1] <- as.integer(counteval)
  cnt[2] <- NA
  names(cnt) <- c("function", "gradient")

  ## Currently only the 'best' solution is returned, but since the
  ## CMA-ES is not a 'global' algorithm in the sense SANN is, all
  ## other solutions are supposed to be 'near' to the best on in the
  ## parameter space.
  res <- list(par=arx[, arindex[1]],
              value=arfitness[1],
              counts=cnt,
              convergence=ifelse(counteval >= stopeval, 1L, 0L),
              message=NULL)
  return(res)
}
