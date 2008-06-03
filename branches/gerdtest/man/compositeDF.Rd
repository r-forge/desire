\name{compositeDF}
\alias{compositeDF}
\alias{compositeDF.call}
\alias{compositeDF.function}
\alias{compositeDF.lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Construct composite desirability functions}
\description{
  Combines a desirability function with an inner function and returns a
  new function object which calculates \eqn{d(f(x))}.
}
\usage{
compositeDF(expr, d, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr}{any valid inner function}
  \item{d}{desirability function}
  \item{\dots}{\dots}
}
\details{
  Currently specialized compositioning functions exists for expressions,
  functions and objects with class \code{lm}.
}
\value{
  A composite desirability function.
}
\author{
  Heike Trautmann \email{trautmann@statistik.uni-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.uni-dortmund.de}
}
\seealso{\code{\link{minimumDI}}}
\examples{
h <- harrington2(-1, 1, 1)
## Calculate h(x^2):
ch1 <- compositeDF(x^2, h)
## Calculate h(f(x))
f <- function(x) 2*x + 3
ch2 <- compositeDF(f, h)
}
\keyword{optimize}
\keyword{multivariate}
