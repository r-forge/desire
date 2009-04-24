\name{plot.desire.function}
\alias{plot.desire.function}
\title{Plot curve of desirability function}
\description{
  TBD
}
\usage{
\method{plot}{desire.function}(x, n = 600, xlim = NULL, ylim = c(0, 1),
              xlab = "Value", ylab = "Desirability", ..., main)
}
\arguments{
  \item{x}{desirability function}
  \item{n}{the number of x values at which to evaluate.}
  \item{xlim}{numeric of length 2; if specified, it serves as the default
    for the range of x values.}
  \item{ylim}{numeric of length 2; if specified, it serves as the default
    for the range of y values.}
  \item{xlab}{x axes label}
  \item{ylab}{y axes label}
  \item{\dots}{arguments passed to \code{lines}.}
  \item{main}{main title of plot}
}
\author{
  Heike Trautmann \email{trautmann@statistik.tu-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.tu-dortmund.de}
}
\keyword{hplot}
