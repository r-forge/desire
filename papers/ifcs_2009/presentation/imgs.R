require(desiRe)

w <- 6.6/2.54
h <- w*3/4

## Harrington
d1 <- harrington1(-1, .1, 0, .5)

pdf(file="h1.pdf", width=w, height=h)
opar <- par(mar=c(2, 2, .5, .5), cex=0.8, lwd=1.5)
plot(d1, xlim=c(-1.5, 2.5), main="")
abline(h=c(0, 1), col=grey(.5), lty=2)
points(c(-1, 0), c(.1, .5), pch=20)
par(opar)
dev.off()

d2a <- harrington2(-1, 1, 1)
d2b <- harrington2(-1, 1, .33)
d2c <- harrington2(-1, 1, 2.5)

pdf(file="h2.pdf", width=w, height=h)
opar <- par(mar=c(2, 2, .5, .5), cex=0.8, lwd=1.5)
plot(d2a, main="")
abline(h=c(0, 1), col=grey(.5), lty=2)
abline(v=c(-1, 1), col=grey(.5), lty=2)
curve(d2b(x), -4, 4, add=TRUE, , n=205, col="blue")
curve(d2c(x), -4, 4, add=TRUE, , n=205, col="red")
par(opar)
dev.off()

## Derringer Suich
d1a <- derringerSuich(c(-Inf, 0, 1, 1, 1))
d1b <- derringerSuich(c(-Inf, 0, 1, 1, .33))
d1c <- derringerSuich(c(-Inf, 0, 1, 1, 2.5))
pdf(file="d1.pdf", width=w, height=h)
opar <- par(mar=c(2, 2, .5, .5), cex=0.8, lwd=1.5)
plot(d1a, main="", , xlim=c(-.25, 1.25))
abline(h=c(0, 1), col=grey(.5), lty=2)
curve(d1b, add=TRUE, n=205, col="blue")
curve(d1c, add=TRUE, n=205, col="red")
par(opar)
dev.off()

d2a <- derringerSuich(c(-1, 0, 2, 1, 1))
d2b <- derringerSuich(c(-1, 0, 2, .33, 2.5))
pdf(file="d2.pdf", width=w, height=h)
opar <- par(mar=c(2, 2, .5, .5), cex=0.8, lwd=1.5)
plot(d2a, main="", , xlim=c(-1.25, 2.25))
abline(h=c(0, 1), col=grey(.5), lty=2)
curve(d2b, add=TRUE, n=205, col="blue")
par(opar)
dev.off()
