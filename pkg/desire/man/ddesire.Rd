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
  Heike Trautmann \email{trautmann@statistik.tu-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.tu-dortmund.de}
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
\examples{
data(Chocolate)

## Fit linear model to data:
m.d90 <- lm(d90 ~ rt + as + I(rt^2) + I(as^2) + rt:as, Chocolate)
m.Fe <- lm(Fe ~ rt + as + I(rt^2) + I(as^2) + rt:as, Chocolate)

## Define desirability functions:
d.d90 <- harrington2(21, 22, 1)
d.Fe <- harrington1(22, 0.8, 28, 0.2)

## Plot density of desirability in rt=30, as=50:
df <- data.frame(rt=30, as=50)
y.Fe <- predict(m.Fe, df)
sigma.Fe <- summary(m.Fe)$sigma

y.d90 <- predict(m.d90, df)
sigma.d90 <- summary(m.d90)$sigma

## Plot curve of density function:
opar <- par(mfrow=c(2,1))
curve(ddesire(x, d.d90, y.d90, sigma.d90), 0, 1, main="d.90", n=202)
curve(ddesire(x, d.Fe, y.Fe, sigma.Fe), 0, 1, main="Fe", n=202)
par(opar)

## Integrate:
integrate(function(x) ddesire(x, d.d90, y.d90, sigma.d90), 0, 1)
integrate(function(x) ddesire(x, d.Fe, y.Fe, sigma.Fe), 0, 1)
}
\keyword{distribution}
