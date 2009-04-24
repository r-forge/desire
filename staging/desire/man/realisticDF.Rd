\name{realisticDF}
\alias{realisticDF}
\alias{realisticDF.desire.function}
\title{Realistic desirability functions}
\description{
  Convert a desirability into a realistic desirability.
}
\usage{
realisticDF(f, ...)
\method{realisticDF}{desire.function}(f, ...)
}
\arguments{
  \item{f}{desirability function}
  \item{\dots}{\dots}
}
\details{
  To construct a composite realistic desirability, first create a
  realistic desirability and then compose it. Doing the opposite is
  currently unsupported. This allows the composition to possibly pass a
  standard deviation which can be deduced from the inner function
  (eg. if the inner function is an object of class \code{lm}).
}
\value{
  A function with the same arguments as \code{x} and \code{sd}, which
  returns the realistic desirability.
}
\author{
  Heike Trautmann \email{trautmann@statistik.tu-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.tu-dortmund.de}
}
\seealso{\code{\link{minimumDI}}}
\keyword{optimize}
\keyword{multivariate}
