require(desire)

w <- 5.2/2.54
h <- w*3/4

## Harrington
d1 <- harrington1(-1, .1, 0, .5)

pdf(file="h1.pdf", width=w, height=h)
opar <- par(mar=c(2, 2, .5, .5), cex=.6)
plot(d1, xlim=c(-1.5, 2.5), main="",  lwd=.5)
points(c(-1, 0), c(.1, .5), pch=20)
par(opar)
dev.off()

d2a <- harrington2(-1, 1, 1)
d2b <- harrington2(-1, 1, .33)
d2c <- harrington2(-1, 1, 2.5)

pdf(file="h2.pdf", width=w, height=h)
opar <- par(mar=c(2, 2, .5, .5), cex=.6)
plot(d2a, main="", lwd=.5)
curve(d2b(x), -4, 4, add=TRUE, lwd=.5, n=205, col="blue")
curve(d2c(x), -4, 4, add=TRUE, lwd=.5, n=205, col="red")
par(opar)
dev.off()

## Derringer Suich
d1a <- derringerSuich(c(-Inf, 0, 1, 1, 1))
d1b <- derringerSuich(c(-Inf, 0, 1, 1, .33))
d1c <- derringerSuich(c(-Inf, 0, 1, 1, 2.5))
pdf(file="d1.pdf", width=w, height=h)
opar <- par(mar=c(2, 2, .5, .5), cex=.6)
plot(d1a, main="", lwd=.5)
curve(d1b, add=TRUE, lwd=.5, n=205, col="blue")
curve(d1c, add=TRUE, lwd=.5, n=205, col="red")
par(opar)
dev.off()

d2a <- derringerSuich(c(-1, 0, 2, 1, 1))
d2b <- derringerSuich(c(-1, 0, 2, .33, 2.5))
pdf(file="d2.pdf", width=w, height=h)
opar <- par(mar=c(2, 2, .5, .5), cex=.6)
plot(d2a, main="", lwd=.5)
curve(d2b, add=TRUE, lwd=.5, n=205, col="blue")
par(opar)
dev.off()
