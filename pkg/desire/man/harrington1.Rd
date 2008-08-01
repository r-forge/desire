\name{harrington1}
\alias{harrington1}
\alias{ddesire.harrington1}
\alias{pdesire.harrington1}
\alias{qdesire.harrington1}
\alias{edesire.harrington1}
\alias{vdesire.harrington1}
\alias{dharrington1}
\alias{pharrington1}
\alias{qharrington1}
\alias{rharrington1}
\alias{eharrington1}
\alias{vharrington1}
\title{One sided Harrington type desirability function}
\description{
  Returns a one sided desirability function of the Harrington type.
}
\usage{
harrington1(y1, d1, y2, d2)
\method{ddesire}{harrington1}(x, f, mean, sd)
\method{pdesire}{harrington1}(q, f, mean, sd)
\method{qdesire}{harrington1}(p, f, mean, sd)
\method{edesire}{harrington1}(f, mean, sd)
\method{vdesire}{harrington1}(f, mean, sd)
dharrington1(x, y1, d1, y2, d2, mean, sd)
pharrington1(q, y1, d1, y2, d2, mean, sd)
qharrington1(p, y1, d1, y2, d2, mean, sd)
rharrington1(n, y1, d1, y2, d2, mean, sd)
eharrington1(y1, d1, y2, d2, mean, sd)
vharrington1(y1, d1, y2, d2, mean, sd)
}
\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilies.}
  \item{n}{number of observations.}
  \item{f}{one sided Harrington type desirability function.}
  \item{y1}{y1}
  \item{d1}{d1}
  \item{y2}{y2}
  \item{d2}{d2}
  \item{mean}{vector of means.}
  \item{sd}{vector of standard deviations.}
}
\value{
  \code{harrington1} returns a function.
}
\references{
  J. Harrington (1965),
  The desirability function.
  \emph{Industrial Quality Control},
  \bold{21}:494-498.
}
\author{
  Heike Trautmann \email{trautmann@statistik.uni-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.uni-dortmund.de}
}
\seealso{\code{\link{harrington2}} for two sided Harrington type
  desirabilities}
\examples{
h <- harrington1(-2, .1, 2, .9)
## Simple checks:
stopifnot(abs(environment(h)$b0 - 0.708167441032245) < 1e-10)
stopifnot(abs(environment(h)$b1 - 0.7710999431401) < 1e-10)
stopifnot(abs(h(-2) -.1) < 1e-10)
stopifnot(abs(h(2) - .9) < 1e-10)

## Plot of desirability:
plot(h)

## d/p/q/r examples:
ddesire(.8, h, 0, 1)
dharrington1(.8, -2, .1, 2, .9, 0, 1)

pdesire(.8, h, 0, 1)
pharrington1(.8, -2, .1, 2, .9, 0, 1)

qdesire(.8, h, 0, 1)
qharrington1(.8, -2, .1, 2, .9, 0, 1)

## Simple convergence plot:
set.seed(42)
x <- rdesire(1e6, h, 0, 1)
est <- cumsum(x)/(1:length(x))
plot(est, type="l")
abline(h=edesire(h, 0, 1), lty=2, col="red")
}
\keyword{optimize}
\keyword{multivariate}
