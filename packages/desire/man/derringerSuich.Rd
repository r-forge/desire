\name{derringerSuich}
\alias{derringerSuich}
\title{Derringer-Suich type desirability function}
\description{
  TBD
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
  If only \code{y} is provided and it is a vector of length 5, a
  \eqn{(l, t, u, \beta_1, \beta_2)} Derringer Suich desirability is
  constructed. Otherwise \code{y}, \code{d} and \code{beta} specify a
  generalized Derringer Suich type desirability.
}
\value{
  \code{derringerSuich} returns a function.
}
\references{TBD}
\author{
  Heike Trautmann \email{trautmann@statistik.uni-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.uni-dortmund.de}
}
\seealso{
  \code{\link{harrington1}} for one sided Harrington type desirabilities and
  \code{\link{harrington2}} for two sided Harrington type desirabilities.1
}
\keyword{optimize}
\keyword{multivariate}
