\name{generationalDistance}
\alias{generationalDistance}
\alias{generalizedSpread}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Quality measures for MCO solutions}
\description{
  Functions to evaulate the quality of the estimated pareto front.
}
\usage{
generationalDistance(x, o)
generalizedSpread(x, o)
}
\arguments{
  \item{x}{Estimated pareto front or an object which has a paretoFront method}
  \item{o}{True pareto front or an object which has a paretoFront method}
}
\details{
  Instead of the pareto front, one can also pass an object for which a
  paretoFront method exists to both methods.
}
\value{
  The respective quality measure.
}
\author{
  Heike Trautmann \email{trautmann@statistik.uni-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.uni-dortmund.de}
}
\examples{
sch1 <- function(x) c(x^2, (x-2)^2)

## Estimate true front:
tf <- nsga2(sch1, 1, 2,
            lower.bounds=0, upper.bounds=1,
            popsize=1000, generations=20)
res <- nsga2(sch1, 1, 2,
             lower.bounds=0, upper.bounds=1,
             popsize=16, generations=c(2, 4, 6, 8, 10, 20, 50))
n <- length(res)
sapply(1:n, function(i) generationalDistance(res[[i]], tf))
sapply(1:n, function(i) generalizedSpread(res[[i]], tf))
}
