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

colramp <- grey(1:200/201)

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
opar <- par(mfrow=c(3, 2))
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
cd.rt     <- compositeDF(fnrt, d.rt)
cd.E      <- compositeDF(m.E, d.E)
cd.d90    <- compositeDF(m.d90, d.d90)
cd.Fe     <- compositeDF(m.Fe, d.Fe)
cd.etaCa  <- compositeDF(m.etaCa, d.etaCa)
cd.tauCa  <- compositeDF(m.tauCa, d.tauCa)

## Realistic Desirabilities:
cd.E      <- compositeDF(m.E, realisticDF(d.E))
cd.d90    <- compositeDF(m.d90, realisticDF(d.d90))
cd.Fe     <- compositeDF(m.Fe, realisticDF(d.Fe))
cd.etaCa  <- compositeDF(m.etaCa, realisticDF(d.etaCa))
cd.tauCa  <- compositeDF(m.tauCa, realisticDF(d.tauCa))

## Combine desirabilities using geometric DI:
di <- geometricDI(cd.rt, cd.E, cd.d90, cd.Fe, cd.etaCa, cd.tauCa,
                  weights=c(10, .1, 5, 1, 1, 1))

## Optimize!
opt <- optim(c(40, 70), di, control=list(fnscale=-1))
opt <- optim(c(40, 70), di, method="BFGS", control=list(fnscale=-1))

opt$par

## Visualize DI:
n <- 101
rt <- seq(20, 45, length.out=n)
as <- seq(5, 90, length.out=n)
sq <- data.frame(rt=rep(rt, each=n), as=rep(as, times=n))

z <- di(sq)
panel <- mkPanel(opt$par[1], opt$par[2], panel.levelplot)
levelplot(z ~ rt * as, sq, panel=panel, cuts=50, col.regions=colramp)

## IndividUalir desirabilities:
res <- NULL
for (dfn in c(cd.rt, cd.E, cd.d90, cd.Fe, cd.etaCa, cd.tauCa)) {
  z <- dfn(sq)
  res <- rbind(res, data.frame(fn=responseName(dfn), rt=sq$rt, as=sq$as, di=z))
}
levelplot(di ~ rt + as | fn, res,
          panel=panel, cuts=50, col.regions=colramp)
