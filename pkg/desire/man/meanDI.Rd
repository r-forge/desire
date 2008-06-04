\name{meanDI}
\alias{meanDI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Weighted mean desirability index}
\description{
  TBD
}
\usage{
meanDI(f, ..., weights = 1)
}
\arguments{
  \item{f,\dots}{\dots}
  \item{weights}{vector of weights. Need not sum to one.}
}
\author{
  Heike Trautmann \email{trautmann@statistik.uni-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.uni-dortmund.de}
}
\value{
  A function object implementing which calculates the weighted mean of
  the desirabilities.
}
\seealso{\code{\link{geometricDI}} and \code{\link{minimumDI}}}
\keyword{optimize}
\keyword{multivariate}
