desireTest <- function() {
  tpar <- par(ask=TRUE)
  on.exit(par(tpar))
  ## Harrington:
  h1 <- harrington2(-5, 2, 2)
  h2 <- harrington1(-2, 0.2, 2, 0.9)
  
  opar <- par(mfcol=c(2,1))
  plot(h1)
  plot(h2)
  par(opar)
  
  ## Derringer Suich:
  ## Taken from Detlefs Dissertation p.40
  f <- derringerSuich(c(-4.8, -3, -1, 0, 1.25, 2.75, 4.6),
                      c(0, 0.2, 0.9, 1, 0.7, 0.2, 0),
                      c(1, 0.5, 4, 0.8, 2, 1))
  
  g <- derringerSuich(c(-5, 0, 5),
                      c(0, 1, 0),
                      c(1, 1))
  
  i <- desirebilityIndex(c(f, g))
  
  plot(f, 800, xlim=c(-6, 6))

  ## Fun with contour()

}
