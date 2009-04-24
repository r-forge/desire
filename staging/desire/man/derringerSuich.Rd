\name{derringerSuich}
\alias{derringerSuich}
\title{Derringer-Suich type desirability function}
\description{
  Returns a function evaluating a (generalized) Derringer-Suich desirability. 
}
\usage{
derringerSuich(y, d, beta)
}
\arguments{
  \item{y}{y}
  \item{d}{d}
  \item{beta}{beta}
}
\details{
  If only \code{y} is provided and it is a vector of length 5, a normal
  \eqn{(l, t, u, \beta_1, \beta_2)} Derringer Suich desirability is
  constructed. Otherwise \code{y}, \code{d} and \code{beta} specify a
  generalized Derringer Suich type desirability.
}
\value{
  \code{derringerSuich} returns a function.
}
\references{
G. Derringer and R. Suich (1980),
Simultaneous optimization of several response variables,
\emph{Journal of Quality Technology},
\bold{Vol. 12, No. 4}:214-219.

D. Steuer (2005),
Statistische Eigenschaften der Multikriteriellen Optimierung mittels Wuenschbarkeiten,
PhD Thesis,
\url{http://hdl.handle.net/2003/20171}.
}
\author{
  Heike Trautmann \email{trautmann@statistik.tu-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.tu-dortmund.de}
}
\seealso{
  \code{\link{harrington1}} for one sided Harrington type desirabilities and
  \code{\link{harrington2}} for two sided Harrington type desirabilities.
}
\examples{
## Simple (l, t, u, b0, b1) Derringer-Suich desirabilities:
d1 <- derringerSuich(c(0, 1, 2, 1, 1))
d2 <- derringerSuich(c(0, 1, 2, 2, 2))
d3 <- derringerSuich(c(0, 1, 2, .5, .5))

## Comparison of their shape:
opar <- par(mfrow=c(3, 1))
plot(d1)
plot(d2)
plot(d3)
par(opar)

## d/p/q/r examples:
ddesire(.2, d1, 0, 1)
pdesire(.5, d1, 0, 1)
qdesire(.8, d1, 0, 1)
}
\keyword{optimize}
\keyword{multivariate}
