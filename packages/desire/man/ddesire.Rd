\name{ddesire}
\alias{ddesire}
\alias{ddesire.default}
\alias{pdesire}
\alias{pdesire.default}
\alias{qdesire}
\alias{qdesire.default}
\alias{rdesire}
\alias{rdesire.default}
\alias{edesire}
\alias{edesire.default}
\alias{vdesire}
\alias{vdesire.default}
\title{Generic Distribution functions for desirabilities}
\description{
  Generic density, distribution, quantile and random number generation functions
  for desirability functions.
}
\usage{
\method{ddesire}{default}(x, f, mean = 0, sd = 1)
\method{pdesire}{default}(q, f, mean = 0, sd = 1)
\method{qdesire}{default}(p, f, mean = 0, sd = 1)
\method{rdesire}{default}(n, f, mean = 0, sd = 1)
\method{edesire}{default}(f, mean, sd)
\method{vdesire}{default}(f, mean, sd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,q}{Vector of quantiles.}
  \item{p}{vector of probabilies.}
  \item{n}{number of observations.}
  \item{f}{desirability function}
  \item{mean}{vector of means.}
  \item{sd}{vector of standard deviations.}
}
\value{
  'ddesire' gives the density, 'pdesire' gives the distribution
  function, 'qdesire' gives the quantile function, and 'rdesire'
  generates random deviates.

  'edesire' and 'vdesire' return the expectation and variance of the
  function.
}
\author{
  Heike Trautmann \email{trautmann@statistik.uni-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.uni-dortmund.de}
}
\note{
  The default implementations for \code{pdesire}, \code{qdesire},
  \code{edesire} and \code{vdesire} are only approximations obtained by
  estimating the desired property from a random sample.
}
\seealso{
  For desirability functions:
  \code{\link{harrington1}} and 
  \code{\link{harrington2}}}
\keyword{distribution}
