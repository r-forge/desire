\name{genShiftedRosenbrock}
\alias{genShiftedRosenbrock}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Single objective test functions}
\description{
  Test and benchmark function generators for single objective
  optimization procedures.
}
\usage{
genShiftedRosenbrock(o, bias)
}
\arguments{
  \item{o}{Location of global optimum.}
  \item{bias}{Function bias.}
}
\value{
  For \code{genShiftedRosenbrock} a function which evaluates to the
  specified shifted Rosenbrock function.
}
\author{
  Olaf Mersmann \email{olafm@statistik.tu-dortmund.de}
}
\seealso{See Also as \code{\link{cma.es}} for a usage example.}
\keyword{optimize}
