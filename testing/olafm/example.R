require(desire)

data(mtcars)

## Derringer-Suich oder Harrington DFs?
use.DS <- FALSE

## Realistische W-Keiten benutzen?
use.RD <- FALSE

if (use.DS) {
  ## Derringer-Suich:
  d.mpg  <- derringerSuich(y=c(-Inf, 15, 20), d=c(1, 1, 0), beta=c(1, 2))
  ## Bsp. für (l, t, u, b1, b2) DS W-keit:
  d.hp   <- derringerSuich(c(80, 160, 400, 1, 2))
  d.qsec <- derringerSuich(y=c(-Inf, 16, 22), d=c(1, 1, 0), beta=c(1, 2))
} else {
  ## Harrington:
  d.mpg  <- harrington1(15, .9, 20, .1)
  d.hp   <- harrington2(140, 180, 1)
  d.qsec <- harrington1(16, .9, 22, .1)
}

## Zugehörige realistische Wünschbarkeit:
rd.mpg  <- realisticDF(d.mpg)
rd.hp   <- realisticDF(d.hp)
rd.qsec <- realisticDF(d.qsec)

## Plot der Wünschbarkeiten
opar <- par(mfrow=c(2,2))
plot(d.mpg, xlab="Miles per Gallon")
plot(d.hp, xlab="Horespower")
plot(d.qsec, xlab="Quarter Mile Time")
par(opar)

gi <- geometricDI(d.mpg, d.hp, d.qsec, weights=c(2, 1, 1))
mi <- minimumDI(d.mpg, d.hp, d.qsec)

Y <- as.matrix(mtcars[,c("mpg", "hp", "qsec")])

## Für den Maserati Bora die beiden DI berechnen:
maserati <- Y["Maserati Bora",]
gi(maserati)
mi(maserati)

## Für alle Fahrzeuge die Wünschbarkeit berechnen und sortiert ausgeben:
sort(apply(Y, 1, gi))
sort(apply(Y, 1, mi))

## Regressionsmodelle für die Zielvariablen:
m.mpg  <- lm(mpg  ~ cyl + disp + drat + wt, mtcars)
m.hp   <- lm(hp   ~ cyl + disp + drat + wt, mtcars)
m.qsec <- lm(qsec ~ cyl + disp + drat + wt, mtcars)

## Zusammengesetze W-keit
if (use.RD) {
  cd.mpg  <- compositeDF(m.mpg , rd.mpg)
  cd.hp   <- compositeDF(m.hp  , rd.hp)
  cd.qsec <- compositeDF(m.qsec, rd.qsec)
} else {
  cd.mpg  <- compositeDF(m.mpg , d.mpg)
  cd.hp   <- compositeDF(m.hp  , d.hp)
  cd.qsec <- compositeDF(m.qsec, d.qsec)
}

## DI über zusammengesetze W-keit
cgi <- geometricDI(cd.mpg, cd.hp, cd.qsec, weights=c(2, 1, 1))
cmi <- minimumDI(cd.mpg, cd.hp, cd.qsec)

## Vergleich DI aus lm und anhand der realen Daten:
cbind(apply(Y, 1, gi), cgi(mtcars))

## Optimiere den Index:
##
## Wichtig ist 'fnscale=-1' damit maximiert wird
##
## Ausserdem nutze ich BFGS, was eigentlich einen gradienten
## voraussetzt, es ist aber das einzige in R implementierte Verfahren,
## dass mit Restriktionen umgehen kann. Bei Nelder Mead werden für
## diesen Anwendungsfall völlig unsinnige Parameter gewählt. Allg. ist
## es ein schlecht gewähltes Beispiel, denn cyl müsste eigentlich
## ganzzahlig sein...
##
## 'cgi' kann auch gegen 'cmi' ausgetauscht werden, wenn man den
## Minimums-Index optimieren möchte. Dann ist die Anwendung von
## L-BFGS-B aber endgültig fragwürdig.
o <- optim(par=c(8, 300, 3.5, 3.5), cgi,
           method="L-BFGS-B",
           lower=c(4, 80, 1, 1),
           upper=c(12, 360, 6, 4),
           control=list(fnscale=-1, trace=1))

print(o)
cgi(o$par)
