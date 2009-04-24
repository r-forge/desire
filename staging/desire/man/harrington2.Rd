\name{harrington2}
\alias{harrington2}
\alias{ddesire.harrington2}
\alias{pdesire.harrington2}
\alias{qdesire.harrington2}
\alias{edesire.harrington2}
\alias{vdesire.harrington2}
\alias{dharrington2}
\alias{pharrington2}
\alias{qharrington2}
\alias{rharrington2}
\alias{eharrington2}
\alias{vharrington2}
\title{Two-sided Harrington type desirability function}
\description{  
  Returns a two sided desirability function of the Harrington type.
  Density, distribution function, quantile function and random number
  generation for the distribution of the two-sided Harrington
  desirability function are computed given a normally distributed
  variable \eqn{Y} with expected value equal to \code{mean} and standard
  deviation equal to \code{sd}.
}
\usage{
harrington2(LSL, USL, n)
\method{ddesire}{harrington2}(x, f, mean, sd)
\method{pdesire}{harrington2}(q, f, mean, sd)
\method{qdesire}{harrington2}(p, f, mean, sd)
dharrington2(x, LSL, USL, n, mean, sd)
pharrington2(q, LSL, USL, n, mean, sd)
qharrington2(p, LSL, USL, n, mean, sd)
rharrington2(ns, LSL, USL, n, mean, sd)
eharrington2(LSL, USL, n, mean, sd)
vharrington2(LSL, USL, n, mean, sd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilies.}
  \item{ns}{number of observations.}
  \item{f}{two-sided Harrington type desirability function.}
  \item{LSL}{Lower Specification Limit of \eqn{Y}.}
  \item{USL}{Upper Specification Limit of \eqn{Y}.}
  \item{n}{Kurtosis parameter of desirability function. Values > 1
    result in smoother shapes around the target value \eqn{T =
      (LSL+USL)/2}. Values < 1 already penalize small target
    deviations.} 
  \item{mean}{vector of means.}
  \item{sd}{vector of standard deviations.}
}

\details{
  \code{harrington2(LSL, USL, n)} is the two-sided desirability function
  of Harrington type (Harrington (1965)). It aims at the specification
  of desired values of a variable \eqn{Y} which has to be optimized
  regarding a target value \eqn{T}.  \eqn{Y} is transformed onto a
  unitless scale to the interval [0,1]. \eqn{LSL} and \eqn{USL} are
  associated with a desirability of \eqn{1/e \approx 0.37}{1/e.,
  approx. 0.37}.  \eqn{LSL} and \eqn{USL} have to be chosen
  symmetrically around the target value \eqn{T}.

  The density and distribution functions of Harrington's two-sided
  desirability function \code{d} given a normally distributed variable
  \eqn{Y} with \eqn{E(Y)=} \code{mean} and \eqn{sd(Y)=}\code{sd} can be
  determined analytically, see Trautmann and Weihs (2006).
}


\value{
  \code{harrington2(LSL, USL, n)} returns a function object of the
  two-sided desirability function of the Harrington type (see example
  below).

  \code{ddesire} / \code{dharrington2} give the density, \code{pdesire}
  / \code{pharrington2} give the distribution function, \code{qdesire} /
  \code{qharrington2} give the quantile function, and \code{rdesire} /
  \code{rharrington2} generate random deviates.  \code{edesire} /
  \code{eharrington2} and \code{vdesire} / \code{vharrington2} compute
  the expected value and the variance of the desirability function for a
  normally distributed random variable \eqn{Y} with
  \eqn{E(Y)=}\code{mean} and \eqn{sd(Y)=}\code{sd}.
}
\references{
  J. Harrington (1965): The desirability function.
  \emph{Industrial Quality Control},
  \bold{21}:494-498.
  
  H. Trautmann, C. Weihs (2006):
  On the Distribution of the Desirability Index using Harrington's Desirability Function.
  \emph{Metrika} \bold{63}(2): 207-213.
}
\author{
  Heike Trautmann \email{trautmann@statistik.tu-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.tu-dortmund.de}
}
\seealso{\code{\link{harrington1}} for one sided Harrington type desirabilities}

\examples{##Assigning the function object to h: 
h <- harrington2(3,7,1) 

## Plot of desirability function: 
plot(h)

## Desirability function of a vector: 
h(seq(2,8,0.1))

## d/p/q/r/e/v examples: 
ddesire(4, h, 0, 1)
dharrington2(4, 3, 7, 1, 0, 1) 

ddesire(4, h, c(0,0.5),c(1,1.5))

pdesire(4, h, 0, 1)
pharrington2(4, 3, 7, 1, 0, 1)

qdesire(0.8, h, 0, 1)
qharrington2(0.8, 3, 7, 1, 0, 1)

rdesire(1e6, h, 0, 1)
rharrington2(1e6, 3, 7, 1, 0, 1)

edesire(h,3,0.5)

vdesire(h,3,0.5)
}
\keyword{optimize}
\keyword{multivariate}
