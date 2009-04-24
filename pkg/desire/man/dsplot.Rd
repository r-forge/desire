\name{dsplot}
\alias{dsplot}
\title{Desirability Plot}
\description{
  TBD
}
\usage{
dsplot(expr, f, from = NULL, to = NULL, n = 101, show.zero = TRUE, interest = NULL, main = "Desirability Plot", sub = NULL, ...)
}
\arguments{
  \item{expr}{an expression written as a function of 'x', or alternatively
          the name of a function which will be plotted.}
  \item{f}{desirability function}
  \item{from,to}{the range over which the function will be plotted.}
  \item{n}{integer; the number of x values at which to evaluate.}
  \item{show.zero}{add dotted line to visualize the origin}
  \item{interest}{vector of interesting points}
  \item{main}{an overall title for the plot}
  \item{sub}{a subtitle for the plot}
  \item{\dots}{parameters passed to low level plot functions}
}
\references{}
\author{
  Heike Trautmann \email{trautmann@statistik.tu-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.tu-dortmund.de}
}
\examples{
d1 <- harrington2(2, 4, 2)
dsplot(x*log(x), d1,
       from=0, to=6, n=507,
       interest=c(2, 3, 4))
}
\keyword{hplot}
