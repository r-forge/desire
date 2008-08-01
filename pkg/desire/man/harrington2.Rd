\name{harrington2}
\alias{harrington2}
\alias{ddesire.harrington2}
\alias{pdesire.harrington2}
\alias{qdesire.harrington2}
\alias{dharrington2}
\alias{pharrington2}
\alias{qharrington2}
\alias{rharrington2}
\alias{eharrington2}
\alias{vharrington2}
\title{Two sided Harrington type desirability functions}
\description{
  Returns a two sided desirability function of the Harrington type.
}
\usage{
harrington2(LSL, USL, n)
\method{ddesire}{harrington2}(x, f, mean, sd)
\method{pdesire}{harrington2}(q, f, mean, sd)
\method{qdesire}{harrington2}(p, f, mean, sd)
dharrington2(x, LSL, USL, n, mean, sd)
pharrington2(q, LSL, USL, n, mean, sd)
qharrington2(p, LSL, USL, n, mean, sd)
rharrington2(ns, LSL, USL, n, mean, sd)
eharrington2(LSL, USL, n, mean, sd)
vharrington2(LSL, USL, n, mean, sd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilies.}
  \item{ns}{number of observations.}
  \item{f}{two sided Harrington type desirability function.}
  \item{LSL}{LSL}
  \item{USL}{USL}
  \item{n}{n}
  \item{mean}{vector of means.}
  \item{sd}{vector of standard deviations.}
}
\value{
  \code{harrington2} returns a function.
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
\seealso{\code{\link{harrington1}} for one sided Harrington type
  desirabilities}
\examples{
h <- harrington2(2, 4, 1)
## Simple checks:
stopifnot(environment(h)$d == 2)
stopifnot(environment(h)$s == 6)
stopifnot(abs(h(2) - h(4)) < 1e-10)

## Plot of desirability:
plot(h)

## d/p/q/r examples:
ddesire(.8, h, 0, 1)
dharrington2(.8, 2, 4, 1, 0, 1)

pdesire(.8, h, 0, 1)
pharrington2(.8, 2, 4, 1, 0, 1)

qdesire(.2, h, 0, 1)
qharrington2(.2, 2, 4, 1, 0, 1)
}
\keyword{optimize}
\keyword{multivariate}
