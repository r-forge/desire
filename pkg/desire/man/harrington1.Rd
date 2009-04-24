\name{harrington1} 
\alias{harrington1}
\alias{ddesire.harrington1}
\alias{pdesire.harrington1} 
\alias{qdesire.harrington1}
\alias{edesire.harrington1}
\alias{vdesire.harrington1}
\alias{dharrington1} 
\alias{pharrington1}
\alias{qharrington1}
\alias{rharrington1}
\alias{eharrington1}
\alias{vharrington1} 
\title{One-sided Harrington type desirability function} 
\description{
  Returns a one-sided desirability function of the Harrington type.
  Density, distribution function, quantile function and random number
  generation for the distribution of the one-sided Harrington
  desirability function are computed given a normally distributed
  variable \eqn{Y} with expected value equal to \code{mean} and standard
  deviation equal to \code{sd}.
}
\usage{
harrington1(y1, d1, y2, d2)
\method{ddesire}{harrington1}(x, f, mean, sd)
\method{pdesire}{harrington1}(q, f, mean, sd)
\method{qdesire}{harrington1}(p, f, mean, sd)
\method{edesire}{harrington1}(f, mean, sd)
\method{vdesire}{harrington1}(f, mean, sd)
dharrington1(x, y1, d1, y2, d2, mean, sd)
pharrington1(q, y1, d1, y2, d2, mean, sd)
qharrington1(p, y1, d1, y2, d2, mean, sd)
rharrington1(n, y1, d1, y2, d2, mean, sd)
eharrington1(y1, d1, y2, d2, mean, sd)
vharrington1(y1, d1, y2, d2, mean, sd)
}
\arguments{
  \item{x,q}{vector of quantiles.}
  \item{p}{vector of probabilies.}
  \item{n}{number of observations.}
  \item{f}{one-sided Harrington type desirability function.}
  \item{y1,d1,y2,d2}{Two values \code{y1} and \code{y2} of variable
    \eqn{Y} with respective desirability values \code{d1} and \code{d2}
    determine the shape of the desirability function.} 
  \item{mean}{vector of expected values of normal distributions.}
  \item{sd}{vector of standard deviations of normal distributions.}
}

\details{
  \code{harrington1(y1, d1, y2, d2)} is the one-sided desirability
  function of Harrington type (Harrington (1965)). It aims at the
  specification of desired values of a variable \eqn{Y} which has to be
  minimized or maximized. \eqn{Y} is transformed onto a unitless scale
  to the interval [0,1].

  Harrington's one-sided desirability function \code{d} given a normally
  distributed variable \eqn{Y} with \eqn{E(Y)=} \code{mean} and
  \eqn{sd(Y)=}\code{sd} has the Double Lognormal Distribution (Holland
  and Ahsanullah (1989)).
}

\value{
  \code{harrington1(y1, d1, y2, d2)} returns a function object of the
  one-sided desirability function of the Harrington type (see example
  below). Values \eqn{b_0} and \eqn{b_1} of the desirability function
  formula are determined.
  
 \code{ddesire} /\code{dharrington1}  give the density, \code{pdesire} /
 \code{pharrington1} give the distribution function, \code{qdesire} /
 \code{qharrington1} give the quantile function, and \code{rdesire} /
 \code{rharrington1} generate random deviates. 
 \code{edesire} / \code{eharrington1} and \code{vdesire} /
 \code{vharrington1} compute the expected value and the variance of the
 desirability function for a normally distributed random variable
 \eqn{Y} with \eqn{E(Y)=}\code{mean} and \eqn{sd(Y)=}\code{sd}. 
}
\references{
  J. Harrington (1965): The desirability function.
  \emph{Industrial Quality Control},
  \bold{21}: 494-498.
  
  B. Holland, M. Ahsanullah (1989): Further Results on the Distribution of Meinhold and Singpurwalla. 
  \emph{The American Statistician} \bold{43} (4): 216-219.
  
  H. Trautmann, C. Weihs (2006):
  On the Distribution of the Desirability Index using Harrington's Desirability Function.
  \emph{Metrika} \bold{63}(2): 207-213.
}
\author{
  Heike Trautmann \email{trautmann@statistik.tu-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.tu-dortmund.de}
}
\seealso{\code{\link{harrington2}} for two sided Harrington type desirabilities}

\examples{
##Assigning the function object to h: 
h <- harrington1(-2, .1, 2, .9) 

## Plot of desirability function: 
plot(h)

## Desirability function of a vector: 
h(seq(-2,2,0.1))

## d/p/q/r/e/v examples: 
ddesire(.8, h, 0, 1)
 dharrington1(.8, -2, .1, 2, .9, 0, 1) 

ddesire(.8, h, c(0,0.5), c(1,1.5))

pdesire(.8, h, 0, 1)
pharrington1(.8, -2, .1, 2, .9, 0, 1)

qdesire(.8, h, 0, 1)
qharrington1(.8, -2, .1, 2, .9, 0, 1)

rdesire(1e6, h, 0, 1)
rharrington1(1e6, -2, .1, 2, .9, 0, 1)

edesire(h,3,0.5)
eharrington1(-2, .1, 2, .9,3,0.5)

vdesire(h,3,0.5)
vharrington1(-2, .1, 2, .9,3,0.5) 

## b_0 and b_1 values: 
environment(h)$b0 
environment(h)$b1
}

\keyword{optimize} 
\keyword{multivariate}
