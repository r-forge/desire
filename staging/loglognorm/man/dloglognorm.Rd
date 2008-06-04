\name{dloglognorm}
\alias{dloglognorm}
\alias{ploglognorm}
\alias{qloglognorm}
\alias{rloglognorm}
\alias{mloglognorm}
\alias{eloglognorm}
\alias{vloglognorm}
\title{The Double Log Normal Distribution}
\description{
  Density, distribution function, quantile function, random generation
  and expected value function for the double log normal distribution
  with mean equal to 'mean' and standard deviation equal to 'sd'.
}
\usage{
dloglognorm(x, mean = 0, sd = 1)
ploglognorm(q, mean = 0, sd = 1)
qloglognorm(p, mean = 0, sd = 1)
rloglognorm(n, mean = 0, sd = 1)
mloglognorm(moment, mean, sd)
eloglognorm(mean, sd)
vloglognorm(mean, sd)
}
\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilites.}
  \item{n}{number of observations.}
  \item{mean}{vector of means.}
  \item{sd}{vector of standard deviations.}
  \item{moment}{vector of moments}
}
\details{
  If 'mean' or 'sd' are not specified they assume the default values of
  '0' and '1', respectively.
}
\value{
  'dloglognorm' gives the density, 'ploglognorm' gives the distribution
  function, 'qloglognorm' gives the quantile function, 'rloglognorm'
  generates random deviates, 'mloglognorm' returns the rth moment,
  'eloglognorm' gives the expected value of the distirbution and
  vloglognorm the variance.
}
\references{
  B. Holland, M. Ahsanullah (1989): Further Resultson the Distribution
  of Meinhold and Singpurwalla, The American Statistician 43 (4), p. 216-219
}
\author{
  Heike Trautmann \email{trautmann@statistik.uni-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.uni-dortmund.de}
}
\examples{
  x <- seq(0, 1, by=0.05)
  ## Several different shapes of the density:
  par(mfrow=c(3, 1))
  curve(dloglognorm(x, -0.2, 0.2), 0, 1, main="DLN(-0.2, 0.2)")
  curve(dloglognorm(x,  0.2, 1.0), 0, 1, main="DLN(0.2, 2.0)")
  curve(dloglognorm(x,  0.2, 1.8), 0, 1, main="DLN(0.2, 2.0)")

  ## Check precision:
  z <- x - pnorm(qnorm(x, .2, 1.0), .2, 1.0)
  max(z)
}
\keyword{distribution}
