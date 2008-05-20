\name{nsga2}
\alias{nsga2}
\title{NSGA II MOEA}
\description{
  Multicriterion optimization algorithm
}
\usage{
nsga2(fn, idim, odim, ..., constraints = NULL, cdim = 0, lower.bounds = rep(-Inf, idim), upper.bounds = rep(Inf, idim), popsize = 100, generations = 100, cprob = 0.5, cdist = 5, mprob = 0.9, mdist = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fn}{Function}
  \item{idim}{Input dimension}
  \item{odim}{Output dimension}
  \item{\dots}{Arguments passed through to 'fn'}
  \item{constraints}{Constraint function}
  \item{cdim}{Constraint dimension}
  \item{lower.bounds}{Lower bound of input}
  \item{upper.bounds}{Upper bound of input}
  \item{popsize}{Size of population}
  \item{generations}{Number of generations to breed}
  \item{cprob}{Crossing probability}
  \item{cdist}{Crossing distribution index}
  \item{mprob}{Mutation probability}
  \item{mdist}{Mutation distribution index}
}
\details{
  TBD
}
\value{
  TBD
}
\references{}
\author{
  Heike Trautmann \email{trautmann@statistik.uni-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.uni-dortmund.de}
}
\keyword{optimize}
