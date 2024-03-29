 
## Helper function:
responseName <- function(di) {
  ev <- environment(di)
  if ("lm" %in% class(ev$expr)) {
    deparse(formula(environment(di)$expr)[[2]])
  } else {
    "rt"
  }
}

mkPanel <- function(x, y, other.panel) {
  panel <- function(...) {
    other.panel(...)
    panel.abline(h=y, v=x, col="red", lty=2)
  }
  return(panel)
}

## Color ramp used for levelplots. Change to your liking:
## colramp <- grey(1:200/201)
colramp <- rgb(r=seq(1, 0.3921, length.out=100),
               g=seq(1, 0.5882, length.out=100),
               b=seq(1, 0.0000, length.out=100))

###############################################################################
## Demo starts here:
require(desiRe)
require(lattice)

data(Chocolate)
summary(Chocolate)

## Construct desirabilities
d.rt    <- derringerSuich(c(-Inf, 30, 45, 1, 1))
d.E     <- derringerSuich(c(-Inf, 3, 4, 1, 1))
d.d90   <- derringerSuich(c(20, 21, 23, 1, 1))
d.Fe    <- derringerSuich(c(-Inf, 20, 30, 1, 1))
d.etaCa <- derringerSuich(c(1, 1.5, 2, 1, 1))
d.tauCa <- derringerSuich(c(5, 8, 10, 1, 1))

## Use linear models as predictor:
m.E     <- lm(E ~ rt + as + I(rt^2) + I(as^2) + rt:as, Chocolate)
m.d90   <- lm(d90 ~ rt + as + I(rt^2) + I(as^2) + rt:as, Chocolate)
m.Fe    <- lm(Fe ~ rt + as + I(rt^2) + I(as^2) + rt:as, Chocolate)
m.etaCa <- lm(etaCa ~ rt + as + I(rt^2) + I(as^2) + rt:as, Chocolate)
m.tauCa <- lm(tauCa ~ rt + as + I(rt^2) + I(as^2) + rt:as, Chocolate)

## Plot of desirabilities
opar <- par(mfrow=c(2, 3))
plot(d.rt, main="Runtime", xlim=c(25, 50))
curve(realisticDF(d.rt)(x, sd=1), n=203, add=TRUE, col="red")

plot(d.E, main="Energy", xlim=c(2, 5))
curve(realisticDF(d.E)(x, sd=.145), n=203, add=TRUE, col="red")

plot(d.d90, main="Particle size", xlim=c(15, 28))
curve(realisticDF(d.d90)(x, sd=1.73), n=203, add=TRUE, col="red")

plot(d.Fe, main="Iron content", xlim=c(15, 35))
curve(realisticDF(d.Fe)(x, sd=.71), n=203, add=TRUE, col="red")

plot(d.etaCa, main="Plastic viscosity", xlim=c(-1, 4))
curve(realisticDF(d.etaCa)(x, sd=.81), n=203, add=TRUE, col="red")

plot(d.tauCa, main="yield value", xlim=c(3, 12))
curve(realisticDF(d.tauCa)(x, sd=.72), n=203, add=TRUE, col="red")
par(opar)

## Goal:
## Minimize runtime while keeping the rest of the parameters 'in spec'.
fnrt <- function(x) {
  if (is.data.frame(x)) {
    return (x$rt)
  } else if (is.matrix(x)){
    return (x[1,])
  } else {
    return (x[1])
  }
}

## Idealistic Desirabilities:
cid.rt     <- compositeDF(fnrt, d.rt)
cid.E      <- compositeDF(m.E, d.E)
cid.d90    <- compositeDF(m.d90, d.d90)
cid.Fe     <- compositeDF(m.Fe, d.Fe)
cid.etaCa  <- compositeDF(m.etaCa, d.etaCa)
cid.tauCa  <- compositeDF(m.tauCa, d.tauCa)

## Realistic Desirabilities:
crd.rt     <- compositeDF(fnrt, d.rt)
crd.E      <- compositeDF(m.E, realisticDF(d.E))
crd.d90    <- compositeDF(m.d90, realisticDF(d.d90))
crd.Fe     <- compositeDF(m.Fe, realisticDF(d.Fe))
crd.etaCa  <- compositeDF(m.etaCa, realisticDF(d.etaCa))
crd.tauCa  <- compositeDF(m.tauCa, realisticDF(d.tauCa))

## Combine desirabilities using geometric DI:
idi <- geometricDI(cid.rt, cid.E, cid.d90, cid.Fe, cid.etaCa, cid.tauCa,
                  weights=c(10, .1, 5, 1, 1, 1))

rdi <- geometricDI(crd.rt, crd.E, crd.d90, crd.Fe, crd.etaCa, crd.tauCa,
                   weights=c(10, .1, 5, 1, 1, 1))

## Optimize!
##
## Here we use BFGS, but in general one should check, that the function
## does not have any kinks, since BFGS requires a (approximated) gradient.
##
iopt <- optim(c(40, 70), idi, method="BFGS", control=list(fnscale=-1))
ropt <- optim(c(40, 70), rdi, method="BFGS", control=list(fnscale=-1))

message("Idealistic solution: ", paste(iopt$par, collapse=", "))
message("Realistic  solution: ", paste(ropt$par, collapse=", "))

## Visualize DI:
##
## For this, we generate a data grid and evaluate the DI for each cell:
n <- 101
rt <- seq(20, 45, length.out=n)
as <- seq(5, 90, length.out=n)
sq <- data.frame(rt=rep(rt, each=n), as=rep(as, times=n))

## This is a helper panel function, which overlays the respective optimum:
panel <- function(...) {
  panel.levelplot(...)
  panel.points(iopt$par[1], iopt$par[2], col="red", pch=19)
  panel.points(ropt$par[1], ropt$par[2], col="blue", pch=19)
}

d <- data.frame(z=c(idi(sq), rdi(sq)),
                rt=sq$rt,
                as=sq$as,
                class=rep(c("Idealistic", "Realistic"), each=nrow(sq)))
levelplot(z ~ rt * as | class, d, cuts=50, col.regions=colramp, panel=panel)

## IndividUalir desirabilities:
res <- NULL
for (dfn in c(cd.rt, cd.E, cd.d90, cd.Fe, cd.etaCa, cd.tauCa)) {
  z <- dfn(sq)
  res <- rbind(res, data.frame(fn=responseName(dfn), rt=sq$rt, as=sq$as, di=z))
}
levelplot(di ~ rt + as | fn, res,
          panel=panel, cuts=50, col.regions=colramp)
