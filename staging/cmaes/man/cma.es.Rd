\name{cma.es}
\alias{cma.es}
\title{Covariance matrix adapting evolutionary strategy}
\description{
  Global optimization procedure using a covariance matrix adapting
  evolutionary strategy.
}
\usage{
  cma.es(par, fn, ..., control = list())
}
\arguments{
  \item{par}{Initial values for the parameters to be optimized over.}
  \item{fn}{A function to be minimized (or maximized), with first
    argument the vector of parameters over which minimization is to take
    place. It should return a scalar result.}
  \item{\dots}{Further arguments to be passed to \code{fn}.}
  \item{control}{A list of control parameters. See \sQuote{Details}.}
}
\details{
  Note that arguments after \code{\dots} must be matched exactly.

  By default this function performs minimization, but it will maximize
  if \code{control$fnscale} is negative. It tries to be a drop in
  replacement for \code{optim}.

  The \code{control} argument is a list that can supply any of the
  following components:
  \describe{
    \item{\code{fnscale}}{An overall scaling to be applied to the value
      of \code{fn} during optimization. If negative,
      turns the problem into a maximization problem. Optimization is
      performed on \code{fn(par)/fnscale}.}
    \item{\code{maxit}}{The maximum number of iterations. Defaults to
      \code{100*D^2}, where \code{D} is the dimension of the parameter space.}
    \item{\code{stopfitness}}{Stop if function value is smaller than or
      equal to \code{stopfitness}. This is the only way for the CMA-ES
      to "converge".}
    \item{\code{sigma}}{Inital variance estimates. Can be a single
      number or a vector of length \code{D}, where \code{D} is the dimension
      of the parameter space.}
  }
}
\value{
   A list with components:
  \item{par}{The best set of parameters found.}
  \item{value}{The value of \code{fn} corresponding to \code{par}.}
  \item{counts}{A two-element integer vector giving the number of calls
    to \code{fn}. The second element is always zero for call
    compatibility with \code{optim}.}
  \item{convergence}{An integer code. \code{0} indicates successful
    convergence. Error codes are
    \describe{
      \item{\code{1}}{indicates that the iteration limit \code{maxit}
      had been reached.}
    }
  }
  \item{message}{Always set to \code{NULL}, provided for call
  compatibility with \code{optim}.}
}
\source{
  The code is based on the \file{purecmaes.m} by N. Hansen.
}
\references{
Hansen, N. (2006). The CMA Evolution Strategy: A Comparing Review. In
  J.A. Lozano, P. Larranga, I. Inza and E. Bengoetxea (eds.). Towards a
  new evolutionary computation. Advances in estimation of distribution
  algorithms. pp. 75-102, Springer;
}
\author{
  Olaf Mersmann \email{olafm@statistik.tu-dortmund.de}
}
\seealso{See Also \code{\link{optim}} for traditional optimization methods.}
\examples{
## Compare performance of different algorithms on the shifted Rosenbrock function:

## Test dimension
n <- 10

## Random optimum in [-50, 50]^n
opt <- runif(n, -50, 50)
bias <- 0
f <- genShiftedRosenbrock(opt, bias)

## Inital parameter values
start <- runif(n, -100, 100)

res.nm <- optim(start, f, method="Nelder-Mead")
res.gd <- optim(start, f, method="BFGS")
res.cg <- optim(start, f, method="CG")
res.sa <- optim(start, f, method="SANN")
res.es <- cma.es(start, f)
}
\keyword{nonlinear}
\keyword{optimize}
