\name{plot.harrington2}
\alias{plot.harrington2}
\title{Plot of two-sided Harrington type desirability function}
\description{
  Plots a two-sided desirability function of the Harrington type. 
}
\usage{
\method{plot}{harrington2}(x, ...)
}
\arguments{
  \item{x}{two-sided Harrington desirability function.}
  \item{\dots}{additional parameters passed to \code{plot}.}
}

\details{Lower Specification Limit (\eqn{LSL}) and Upper Specification Limit
  (\eqn{USL}) are visualized. The default range of the x-axis is selected automatically.}
  
\value{R-graphics plot object of two-sided Harrington desirability function.}

\references{
  J. Harrington (1965): The desirability function.
  \emph{Industrial Quality Control},
  \bold{21}: 494-498.}


\author{
  Heike Trautmann \email{trautmann@statistik.tu-dortmund.de},
  Detlef Steuer \email{steuer@hsu-hamburg.de} and
  Olaf Mersmann \email{olafm@statistik.tu-dortmund.de}
}

\seealso{\code{\link{harrington2}} for two sided Harrington type desirabilities, 
\code{\link{plot}}}

\examples{##Assigning a function object to h: 
h <- harrington2(3,7,1) 

## Plot of desirability function: 
plot(h)

plot(harrington2(4,8,1.5))}
