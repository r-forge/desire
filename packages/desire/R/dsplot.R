##
## dsplot.R - desirability plot
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <olafm@statistik.uni-dortmund.de>
##

dsplot <- function(expr, f,
                   from=NULL, to=NULL, n=101,
                   interest=NULL,
                   main="Desirability Plot", sub=NULL,...) {
  sexpr <- substitute(expr)
  if (is.name(sexpr)) {
    fcall <- paste(sexpr, "(x)")
    expr <- parse(text = fcall)
  } else {
    if (!(is.call(sexpr) && match("x", all.vars(sexpr), nomatch = 0))) 
      stop("'expr' must be a function or an expression containing 'x'")
    expr <- sexpr
  }

  ## Evaluate expression:
  x <- seq(from, to, length.out=n)
  y <- eval(expr, envir = list(x = x), enclos = parent.frame())

  if (length(interest) > 0) {
    xi <- interest
    yi <- eval(expr, envir=list(x=xi), enclos=parent.frame())
    di <- f(yi)
  }

  ## Save par and restore on exit
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))

  ## Allocate space for main/sub title
  oma <- c(0, 0, 0, 0)
  if (!is.null(main)) oma[3] <- 2
  if (!is.null(sub)) oma[1] <- 2
  par(oma=oma)

  ## Layout of plots: 1/3 desirability, 2/3 expression
  layout(matrix(c(2, 1), ncol=2),
         widths=c(1/3, 2/3),
         heights=c(1,1))

  ## Expression plot
  par(mar=c(5, 3, 1, 1))
  plot(x, y, type="l", ...)

  if (length(interest) > 0) {
    xmin <- par("usr")[1]
    ymin <- par("usr")[3]
    segments(xi, ymin, xi, yi, col="red", lty=2)
    segments(xmin, yi, xi, yi, col="red", lty=2)
  }
  ## Desirability plot
  par(mar=c(5, 2, 1, 0))
  y <- seq(min(y), max(y), length.out=n)
  plot(1-f(y), y,
       xlab=substitute(d(sexpr)), ylab=sexpr,
       type="l", xaxt="n", yaxt="n", ...)
  axis(4, labels=FALSE)
  tck <- axTicks(1)  
  axis(1, at=tck, labels=1-tck)

  if (length(interest) > 0) {
    xmax <- par("usr")[2]
    segments(xmax, yi, 1-di, yi, col="red", lty=2)
  }
  
  ## Global y axis label
  mtext(sexpr, side=2, line=1)

  ## Add main/sub title to plot:
  if (!is.null(main))
    mtext(main, side=3, outer=TRUE)
  if (!is.null(sub))
    mtext(sub, side=1, outer=TRUE)
  invisible()
}
