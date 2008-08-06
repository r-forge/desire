require(desire)
require(mco)

schoko <- read.table("useR_2008-presentation.txt", sep=";", dec=",", head=TRUE)
schoko$RT <- NULL
schoko$srun <- NULL

summary(schoko)

## Use linear models as predictor:

## Energie consumed
m.E     <- lm(E ~ rt + as + I(rt^2) + I(as^2) + rt:as, schoko)
## paricle size
m.d90   <- lm(d90 ~ rt + as + I(rt^2) + I(as^2) + rt:as, schoko)
## iron content
m.Fe    <- lm(Fe ~ rt + as + I(rt^2) + I(as^2) + rt:as, schoko)
## casson plastic viscosity (characterizes melting behavior)
m.etaCa <- lm(etaCa ~ rt + as + I(rt^2) + I(as^2) + rt:as, schoko)
## casson yield value (characterizes melting behavior)
m.tauCa <- lm(tauCa ~ rt + as + I(rt^2) + I(as^2) + rt:as, schoko)

d.rt    <- derringerSuich(c(-Inf, 30, 45, 1, 1))
d.E     <- derringerSuich(c(-Inf, 3, 4, 1, 1))
d.d90   <- derringerSuich(c(20, 21, 23, 1, 1))
d.Fe    <- derringerSuich(c(-Inf, 20, 30, 1, 1))
d.etaCa <- derringerSuich(c(1, 1.5, 2, 1, 1))
d.tauCa <- derringerSuich(c(5, 8, 10, 1, 1))

## Plot of desirabilities
opar <- par(mfrow=c(3, 2))
plot(d.rt, main="Runtime")
plot(d.E, main="Energy")
plot(d.d90, main="Particle size")
plot(d.Fe, main="Iron content")
plot(d.etaCa, main="Plastic viscosity")
plot(d.tauCa, main="yield value")
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

cd.rt     <- compositeDF(fnrt, d.rt)
cd.E      <- compositeDF(m.E, d.E)
cd.d90    <- compositeDF(m.d90, d.d90)
cd.Fe     <- compositeDF(m.Fe, d.Fe)
cd.etaCa  <- compositeDF(m.etaCa, d.etaCa)
cd.tauCa  <- compositeDF(m.tauCa, d.tauCa)

## Combine desirabilities using geometric DI:
di <- geometricDI(cd.rt, cd.E, cd.d90, cd.Fe, cd.etaCa, cd.tauCa,
                  weights=c(10, .1, 5, 1, 1, 1))

## Optimize!
of <- optim(c(40, 70), di, control=list(fnscale=-1))

## Use NSGA-II
fn <- function(x) {
  ##   df <- data.frame(rt=x[1], as=x[2])
  ##   rt <- x[1]
  ##   E <- predict(m.E, df)
  ##   d90 <- predict(m.d90, df)
  ##   Fe <- predict(m.Fe, df)
  ##   etaCa <- predict(m.etaCa, df)
  ##   tauCa <- predict(m.tauCa, df)  
  ##   c(rt, E, (21-d90)^2, Fe, (1.5 - etaCa)^2, (8 - tauCa)^2)
  
  -c(cd.rt(x), cd.E(x), cd.d90(x), cd.Fe(x), cd.etaCa(x), cd.tauCa(x))
}

## res <- nsga2(fn, 2, 6,
##             lower.bounds=c(20, 5), upper.bounds=c(45, 90),
##             popsize=20, generations=100)
## colnames(res$value) <- c("rt", "E", "d90", "Fe", "etaCa", "tauCa")
res <- dget(res)

## Visualize DI:
rt <- seq(20, 45, length.out=100)
as <- seq(5, 90, length.out=100)

z <- outer(rt, as, function(x,y) di(data.frame(rt=x, as=y)))
image(rt, as, z, col=grey(seq(0, 1, length.out=200)))

## Mark DI optimum
abline(v=of$par[1], h=of$par[2], col="red", lty=2)

## Mark NSGA-II population:
points(res$par, col="blue", pch=20)

## Individual desirabilities:
opar <- par(mfrow=c(3, 2))
z <- outer(rt, as, function(x, y) cd.rt(data.frame(rt=x, as=y)))
image(rt, as, z, col=grey(seq(0, 1, length.out=200)), main="RT")
abline(v=of$par[1], h=of$par[2], col="red", lty=2)

z <- outer(rt, as, function(x, y) cd.E(data.frame(rt=x, as=y)))
image(rt, as, z, col=grey(seq(0, 1, length.out=200)), main="E")
abline(v=of$par[1], h=of$par[2], col="red", lty=2)

z <- outer(rt, as, function(x, y) cd.d90(data.frame(rt=x, as=y)))
image(rt, as, z, col=grey(seq(0, 1, length.out=200)), main="d90")
abline(v=of$par[1], h=of$par[2], col="red", lty=2)

z <- outer(rt, as, function(x, y) cd.Fe(data.frame(rt=x, as=y)))
image(rt, as, z, col=grey(seq(0, 1, length.out=200)), main="Fe")
abline(v=of$par[1], h=of$par[2], col="red", lty=2)

z <- outer(rt, as, function(x, y) cd.etaCa(data.frame(rt=x, as=y)))
image(rt, as, z, col=grey(seq(0, 1, length.out=200)), main="etaCa")
abline(v=of$par[1], h=of$par[2], col="red", lty=2)

z <- outer(rt, as, function(x, y) cd.tauCa(data.frame(rt=x, as=y)))
image(rt, as, z, col=grey(seq(0, 1, length.out=200)), main="tauCa")
abline(v=of$par[1], h=of$par[2], col="red", lty=2)
par(opar)

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
# og <- optim(c(40, 70), rdi, control=list(fnscale=-1))
og <- dget("og.txt")

## Deterministic/Idealized desirability of solution is 0!
di(og$par)

# rz <- outer(rt, as, function(x,y) rdi(data.frame(rt=x, as=y)))
rz <- dget("rz.txt")

image(rt, as, rz, col=grey(seq(0, 1, length.out=200)))
abline(v=og$par[1], h=og$par[2], col="red", lty=2)
abline(v=of$par[1], h=of$par[2], col="blue", lty=2)
