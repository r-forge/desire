##
## plot.R - plot functions for desirebility functions
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##


plot.desire.function <- function(f, n=600,
                                 xlim=NULL, ylim=c(0, 1),
                                 xlab="Value", ylab="Desirebility",
                                 ...) {
  if (is.null(xlim))
    xlim <- attr(f, "y.range")
  plot.new()
  plot.window(xlim, ylim)
  box()
  axis(1)
  axis(2)
  title(main=paste("Type: ", attr(f, "desire.type")),
        xlab=xlab, ylab=ylab)

  ## Don't use xlim as range. Use 'real' range of x axis.
  xrng <- par("usr")[1:2]
  x <- seq(xrng[1], xrng[2], length.out=n)
  y <- f(x)
  abline(h=c(0, 1), col="grey", lty=2)
  lines(x, y, ...)
}
