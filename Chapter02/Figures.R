## 2.1
N = 10^4
x = runif(N)
x1 = x[-1]
x2 = x[-N]

par(mfrow=c(1, 3))
hist(x)
plot(x1, x2, pch=20, cex=0.5)
acf(x)


## 2.2
N = 10^4
x = runif(N)
y1 = -log(1-x)
y2 = rexp(N)
M = max(y1, y2)
X = seq(0, M, 0.01)
Y = dexp(X)

par(mfrow=c(1,2))
breaks = seq(0, M+0.2, 0.2)
hist(y1, freq=F, main="Exp from Uniform", breaks=breaks)
lines(X, Y, type='l', col='red')
hist(y2, freq=F, main="Exp from R", breaks=breaks)
lines(X, Y, type='l', col='red')


## 2.3
N = 10^4
n = 6
p = 0.3
y = rgamma(N, n, scale=(1-p)/p)
X_y = rpois(N, y)

X = seq(0, 60, 1)
hist(X_y, freq=FALSE, breaks=seq(0, 60, 2), main="Neg(6, 3) Mixture Representation")
lines(X, dnbinom(X, n, p), lwd=2, col="sienna")


## 2.4 
N = 2500
f = function(x) { dbeta(x, 2.7, 6.3) }
X = seq(0, 1, 0.001)
U = runif(N)

X1 = runif(N)
g1 = function(x) { dunif(x) }
f_g1 <- function(x) { f(x) / g1(x) }
M1 = optimize(f_g1, c(0, 1), maximum=TRUE)$objective
Y1 = M1 * U * g1(X1)
AcceptX1 = X1[Y1 <= f(X1)]
AcceptY1 = Y1[Y1 <= f(X1)]
RejectX1 = X1[Y1 > f(X1)]
RejectY1 = Y1[Y1 > f(X1)]

X2 = rbeta(N, 2, 6)
g2 = function(x) { dbeta(x, 2, 6) }
f_g2 <- function(x) { f(x) / g2(x) }
M2 = optimize(f_g2, c(0, 1), maximum=TRUE)$objective
Y2 = M2 * U * g2(X2)
AcceptX2 = X2[Y2 <= f(X2)]
AcceptY2 = Y2[Y2 <= f(X2)]
RejectX2 = X2[Y2 > f(X2)]
RejectY2 = Y2[Y2 > f(X2)]

par(mfrow=c(1, 2))

plot(RejectX1, RejectY1, xlim=c(0, 1), pch=20, cex=0.5, col="gray")
points(AcceptX1, AcceptY1, pch=20, cex=0.5)
lines(X, ylim=c(0, M1), M1 * g1(X), col="yellowgreen", lwd=2)
lines(X, f(X), col="red", lwd=2)

plot(RejectX2, RejectY2, xlim=c(0, 1), pch=20, cex=0.5, col="gray")
points(AcceptX2, AcceptY2, pch=20, cex=0.5)
lines(X, M2 * g2(X), col="yellowgreen", lwd=2)
lines(X, f(X), type='l', col="red", lwd=2)