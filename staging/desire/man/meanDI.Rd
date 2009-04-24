\name{meanDI}
\alias{meanDI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Weighted Mean Desirability Index}
\description{Computes the weighted mean of a number of desirability functions. 
 }
\usage{
meanDI(f, ..., weights = 1)
}
\arguments{
  \item{f,\dots}{desirability functions.}
  \item{weights}{vector of weights. Weights do not need to sum to one.}
}


\details{The Desirability Index was introduced by Harrington (1965), and the concept was extended by Derringer and Suich (1980).
It is a means for multicriteria (quality) optimization in industrial quality management. All desirability functions of the 
quality criteria are combined into a univariate global quality criterion which has to be optimized. The Weighted Mean Desirability Index is related
to the concept of utility functions.

The function can be used for Harrington as well as Derringer and Suich desirability functions.}


\value{\code{meanDI(f, ..., weights)} returns  a function object of the Weighted Mean Desirability Index.
}

\references{
  J. Harrington (1965): The desirability function.
  \emph{Industrial Quality Control},
  \bold{21}: 494-498.
  
  G.C. Derringer, D. Suich (1980): Simultaneous optimization of several response variables. 
  \emph{Journal of Quality Technology} \bold{12} (4): 214-219.
  
}

\author{
  Heike Trautmann \email{trautmann@statistik.tu-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.tu-dortmund.de}
}

\seealso{\code{\link{harrington1}} and \code{\link{harrington2}} for Harrington type desirability functions;
         \code{\link{derringerSuich}} for desirability functions of Derringer and Suich;
         \code{\link{geometricDI}},\code{\link{minimumDI}} for other types of Desirability indices.}

\examples{
h1 <- harrington1(-2, .9, 2, .1)
h2 <- harrington2(0, 2, 2)

di <- meanDI(h1, h2,weights=c(0.2,0.8))
di(c(0, 1))

## Desirability Index of vector input:
h <- harrington2(3,7,1)
g <- harrington1(-2, .1, 2, .9) 

d <- meanDI(h, g,weights=c(0.3,0.7))

m <- matrix(c(seq(2, 8, 0.1), seq(-2, 4, 0.1)), ncol=2, byrow=FALSE)
apply(m, 1, d) 


%## Plot of the surface of di:
%x1 <- seq(-3, 3, length.out=20)
%x2 <- seq(-1, 3, length.out=40)
%z <- outer(x1, x2, function(x, y) di(cbind(x, y)))
%image(x1, x2, z, col=grey(seq(0, 1, length.out=120)))
}
\keyword{optimize}
\keyword{multivariate}
