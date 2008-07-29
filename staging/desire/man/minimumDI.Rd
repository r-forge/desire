\name{minimumDI}
\alias{minimumDI}
\title{Minimum desirability Index}
\description{
  TBD
}
\usage{
minimumDI(f, ...)
}
\arguments{
  \item{f,\dots}{desirability functions}
}
\author{
  Heike Trautmann \email{trautmann@statistik.uni-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.uni-dortmund.de}
}
\seealso{\code{\link{geometricDI}}}
\examples{
h1 <- harrington1(-2, .9, 2, .1)
h2 <- harrington2(0, 2, 2)

di <- minimumDI(h1, h2)
di(c(0, 1))

## Plot of the surface of di:
x1 <- seq(-3, 3, length.out=20)
x2 <- seq(-1, 3, length.out=40)
z <- outer(x1, x2, function(x, y) di(cbind(x, y)))
image(x1, x2, z, col=grey(seq(0, 1, length.out=120)))
}
\keyword{optimize}
\keyword{multivariate}
