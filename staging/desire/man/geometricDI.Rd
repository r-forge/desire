\name{geometricDI}
\alias{geometricDI}
\title{Geometric Mean Desirability Index}
\description{
  Computes the weighted geometric mean of a number of desirability functions.
}
\usage{
geometricDI(f, ..., weights)
}
\arguments{
  \item{f, \dots}{desirability functions}
  \item{weights}{vector of weights}
}

\details{
  The Desirability Index was introduced by Harrington (1965), and the
  concept was extended by Derringer and Suich (1980).  It is a means for
  multicriteria (quality) optimization in industrial quality
  management. All desirability functions of the quality criteria are
  combined into a univariate global quality criterion in [0,1] which has
  to be optimized.

  The function can be used for Harrington as well as Derringer and Suich
  desirability functions.
}

\value{
  \code{geometricDI(f, ..., weights)} returns a function object of
  the Geometric Mean Desirability Index.
}

\references{
  J. Harrington (1965): The desirability function.
  \emph{Industrial Quality Control},
  \bold{21}: 494-498.
  
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
         \code{\link{derringerSuich}} for desirability functions of Derringer and Suich;
         \code{\link{minimumDI}},\code{\link{meanDI}} for other types of Desirability indices.}
         
\examples{
h1 <- harrington1(-2, .9, 2, .1)
h2 <- harrington2(0, 2, 2)

di <- geometricDI(h1, h2, weights=c(1/3, 2/3))
di(c(0, 1))

## Desirability Index of vector input:
h <- harrington2(3,7,1)
g <- harrington1(-2, .1, 2, .9) 

d <- geometricDI(h, g, weights=c(0.5, 0.5))

m <- matrix(c(seq(2, 8, 0.1), seq(-2, 4, 0.1)), ncol=2, byrow=FALSE)
apply(m, 1, d) 
}

\keyword{optimize}
\keyword{multivariate}
