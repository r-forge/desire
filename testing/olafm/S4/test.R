source("generics.R")
source("desirability.R")
source("harrington1.R")
source("harrington2.R")
source("random.R")
source("stochastic.R")

test.h1 <- function() {
  f <- harrington1(1, .2, 3.25, .49997)
  g <- harrington1(1, .49997, 3.25, .2)

  opar <- par(mfrow=c(2, 1))
  plot(f)
  plot(g)
  par(opar)
}

test.h2 <- function() {
  opar <- par(mfrow=c(2, 1))
  for (n in c(1, 3.5)) {    
    ## Broken density!
    f <- harrington2(3, 9, n)
    plot(f)
    ## curve(ddesirability(x, f, 5, 1), 0, 1, col="black", ylim=c(0, 5))
    ## curve(ddesirability(x, f, 4, 1), 0, 1, col="blue", add=TRUE)
    ## curve(ddesirability(x, f, 6, 4), 0, 1, col="red", add=TRUE)
  }
  par(opar)
}


test.sd <- function() {
  h <- harrington2(3, 7, 3.5)
  s <- sdesire(h, 1)
  curve(h(x), 2, 8, ylim=c(0, 1))
  curve(s(x), 2, 8, add=TRUE, col="red")
}

