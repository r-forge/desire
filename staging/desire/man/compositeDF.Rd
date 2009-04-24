\name{compositeDF}
\alias{compositeDF}
\alias{compositeDF.call}
\alias{compositeDF.function}
\alias{compositeDF.lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Construct composite desirability functions}
\description{
  Combines a desirability function with an inner function and returns a
  new function object which calculates \eqn{d(f(x))}.
}
\usage{
compositeDF(expr, d, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr}{any valid inner function}
  \item{d}{desirability function}
  \item{\dots}{additional arguments passed to function}
}
\details{
  Currently specialized for expressions,
  functions and objects with class \code{lm}.
}
\value{
  A function object of a composite desirability function.
}

\references{
  J. Harrington (1965): The desirability function.
  \emph{Industrial Quality Control}, \bold{21}: 494-498.
  
  G.C. Derringer, D. Suich (1980): Simultaneous optimization of several response variables. 
  \emph{Journal of Quality Technology} \bold{12} (4): 214-219.
  
  D. Steuer (2005):
  Statistische Eigenschaften der Multikriteriellen Optimierung mittels Wuenschbarkeiten.
  \emph{Dissertation}, Dortmund University of Technology, http://hdl.handle.net/2003/20171.
  
  H. Trautmann, C. Weihs (2006):
  On the Distribution of the Desirability Index using Harrington's Desirability Function.
  \emph{Metrika} \bold{63}(2): 207-213.
}


\author{
  Heike Trautmann \email{trautmann@statistik.tu-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.tu-dortmund.de}
}
\seealso{\code{\link{harrington1}} and \code{\link{harrington2}} for Harrington type desirability functions;
         \code{\link{derringerSuich}} for desirability functions of Derringer and Suich;}
\examples{
h <- harrington2(-1, 1, 1)
## Calculate h(x^2):
ch1 <- compositeDF(x^2, h)
ch1(0.5)
ch1(c(0.2,0.5,0.7))

## Calculate h(f(x))
f <- function(x) 2*x + 3
ch2 <- compositeDF(f, h)
ch2(0.3)
ch2(c(0.3,0.35,0.9))
}
\keyword{optimize}
\keyword{multivariate}
