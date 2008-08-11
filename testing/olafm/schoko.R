require(desire)

schoko <- read.table("schoko.txt", sep=";", dec=",", head=TRUE)
schoko$RT <- NULL
schoko$srun <- NULL

summary(schoko)

m.E     <- lm(E ~ rt + as + I(rt^2) + I(as^2) + rt:as, schoko)
m.d90   <- lm(d90 ~ rt + as + I(rt^2) + I(as^2) + rt:as, schoko)
m.Fe    <- lm(Fe ~ rt + as + I(rt^2) + I(as^2) + rt:as, schoko)
m.etaCa <- lm(etaCa ~ rt + as + I(rt^2) + I(as^2) + rt:as, schoko)
m.tauCa <- lm(tauCa ~ rt + as + I(rt^2) + I(as^2) + rt:as, schoko)

d.rt    <- derringerSuich(c(-Inf, 30, 45, 1, 1))
d.E     <- derringerSuich(c(-Inf, 3, 4, 1, 1))
d.d90   <- derringerSuich(c(20, 21, 23, 1, 1))
d.Fe    <- derringerSuich(c(-Inf, 20, 30, 1, 1))
d.etaCa <- derringerSuich(c(1, 1.5, 2, 1, 1))
d.tauCa <- derringerSuich(c(5, 8, 10, 1, 1))

## Plot der im Artikel angegebenen Wünschbarkeiten:
# opar <- par(mfrow=c(3, 2))
# for (d in c(d.rt, d.E, d.d90, d.Fe, d.etaCa, d.tauCa))
#   plot(d)
# par(opar)

fnrt <- function(x) {
  if (is.data.frame(x)) {
    return (x$rt)
  } else if (is.matrix(x)){
    return (x[1,])
  } else {
    return (x[1])
  }
}

cd.rt     <- compositeDF(fnrt, d.rt)
cd.E      <- compositeDF(m.E, d.E)
cd.d90    <- compositeDF(m.d90, d.d90)
cd.Fe     <- compositeDF(m.Fe, d.Fe)
cd.etaCa  <- compositeDF(m.etaCa, d.etaCa)
cd.tauCa  <- compositeDF(m.tauCa, d.tauCa)

di <- geometricDI(cd.rt, cd.E, cd.d90, cd.Fe, cd.etaCa, cd.tauCa,
                  weights=c(10, .1, 5, 1, 1, 1))

f <- function(x) 
  di(data.frame(rt=x[1], as=x[2]))

of <- optim(c(40, 70), di, control=list(fnscale=-1))

rt <- seq(30, 45, length.out=100)
as <- seq(40, 70, length.out=100)

z <- outer(rt, as, function(x,y) di(data.frame(rt=x, as=y)))
image(rt, as, z, col=grey(seq(0, 1, length.out=200)))
abline(v=of$par[1], h=of$par[2], col="red", lty=2)

## Realistische Version (ACHTUNG: Rechnet ewig wegen DS!)
rcd.E      <- compositeDF(m.E, realisticDF(d.E))
rcd.d90    <- compositeDF(m.d90, realisticDF(d.d90))
rcd.Fe     <- compositeDF(m.Fe, realisticDF(d.Fe))
rcd.etaCa  <- compositeDF(m.etaCa, realisticDF(d.etaCa))
rcd.tauCa  <- compositeDF(m.tauCa, realisticDF(d.tauCa))

rdi <- geometricDI(cd.rt, rcd.E, rcd.d90, rcd.Fe, rcd.etaCa, rcd.tauCa,
                   weights=c(10, .1, 5, 1, 1, 1))

## Interessant ist, dass die Lage des Optimums sich deutlich verändert
## wenn man die realistischen W-Keiten optimiert. 
og <- optim(c(40, 70), rdi, control=list(fnscale=-1))

## Das Optimum hat, bezogen auf die 'normalen' W-Keiten sogar eine Wünschbarkeit von 0:
di(og$par)

rz <- outer(rt, as, function(x,y) rdi(data.frame(rt=x, as=y)))
image(rt, as, rz, col=grey(seq(0, 1, length.out=200)))
abline(v=og$par[1], h=og$par[2], col="red", lty=2)
abline(v=of$par[1], h=of$par[2], col="blue", lty=2)
