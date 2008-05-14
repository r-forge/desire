\name{dloglognorm}
\alias{dloglognorm}
\alias{ploglognorm}
\alias{qloglognorm}
\alias{rloglognorm}
\alias{eloglognorm}
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
eloglognorm(mean = 0, sd = 1)
}
\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilites.}
  \item{n}{number of observations.}
  \item{mean}{vector of means.}
  \item{sd}{vector of standard deviations.}
}
\details{
  If 'mean' or 'sd' are not specified they assume the default values of
  '0' and '1', respectively.
}
\value{
  'dloglognorm' gives the density, 'ploglognorm' gives the distribution
  function, 'qloglognorm' gives the quantile function, 'rloglognorm'
  generates random deviates and 'eloglognorm' gives the expected value
  of the distirbution.  
}
\references{ FIXME: Erstes auftauchen? }
\author{
  Heike Trautmann \email{trautmann@statistik.uni-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.uni-dortmund.de}
}
\keyword{distribution}
