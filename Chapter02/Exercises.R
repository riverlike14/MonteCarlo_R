## 2.2
### (a)
supp = seq(-12, 12, 0.01)
N = 10000
U = runif(N)
f = dlogis
Finv = function(x) { -log((1-x)/x) }
X = Finv(U)
Y = rlogis(N)

par(mfrow=c(1,2))
hist(X, freq=FALSE, breaks=seq(-12, 12, 0.5), main="Logistic from Uniform")
lines(supp, f(supp), col="red", lwd=2)
hist(Y, freq=FALSE, breaks=seq(-12, 12, 0.5), main="Logistic from R")
lines(supp, f(supp), col="red", lwd=2)


### (b)
supp = seq(-12, 12, 0.01)
N = 10000
U = runif(N)
f = dcauchy
Finv = function(x) { tan(pi * (x - 1/2)) }
X = Finv(U); X = X[-12 < X & X < 12]
Y = rcauchy(N); Y = Y[-12 < Y & Y < 12]

par(mfrow=c(1,2))
hist(X, freq=FALSE, breaks=seq(-12, 12, 0.5), main="Cauchy from Uniform")
lines(supp, f(supp), col="red", lwd=2)
hist(Y, freq=FALSE, breaks=seq(-12, 12, 0.5), main="Cauchy from R")
lines(supp, f(supp), col="red", lwd=2)

## 2.3
### (b)
N = 10000
test1 = function(N) {
  U = matrix(runif(12 * N), nrow=12) - 0.5
  Z = apply(U, 2, sum)
}
test2 = function(N) {
  U1 = runif(N)
  U2 = runif(N)
  X1 = sqrt(-2 * log(U1)) * cos(2 * pi * U2)
}
par(mfrow=c(1, 2))
hist(test1(N), freq=FALSE, breaks=seq(-6, 6, 0.2), main="CLT-normal generator")
hist(test2(N), freq=FALSE, breaks=seq(-6, 6, 0.2), main="Box-Muller algorithm")

system.time(test1(N)); system.time(test2(N))

### (c)
test3 = function(N) {
  rnorm(N)
}
hist(test3(N), freq=FALSE)
system.time(test3(N))


## 2.4
### (c)
N = 10000
Sigma = cov(matrix(rnorm(30), nrow=10))

test1 = function(N) {
  A = t(chol(Sigma))
  X = matrix(runif(3 * N), nrow=N)
  X = X %*% A
}
test2 = function(N) {
  library(mnormt)
  rmnorm(N, varcov=Sigma)
}

system.time(test1(N)); system.time(test2(N))



## 2.12
### (a)
a = 3
b = 5
N = 10^4

test1 = function(N, a, b) {
  U = matrix(runif(a * N), nrow=a)
  X = b * apply(-log(U), 2, sum)
}
test2 = function(N, a, b) {
  X = rgamma(N, a, scale=b)
}
system.time(test1(N, a, b)); system.time(test2(N, a, b))

test1 = function(N, a, b) {
  U = matrix(runif((a+b) * N), nrow=a+b)
  X = -log(U)
  X = apply(X[1:a,], 2, sum) / apply(X, 2, sum)
}
test2 = function(N, a, b) {
  rbeta(N, a, b)
}
system.time(test1(N, a, b)); system.time(test2(N, a, b))


## 2.14
### (b)
N = 1000
r = 10

par(mfrow=c(1,3))
for (p in c(0.01, 0.1, 0.5)) {
  tmean = r*(1-p)/p
  tstd = sqrt(r*(1-p)/p^2)
  mrange = round(seq(max(0, tmean - 3*tstd), tmean + 3*tstd, 1))
  prob = pnbinom(mrange, r, p)
  X = rep(0, N)
  
  for (i in 1:N) {
    u = runif(1)
    X[i] = mrange[1] + sum(prob < u)
  }
  
  hist(X, freq=FALSE, main=paste("Negative Binomial with p=", p, sep=""))
  lines(dnbinom(seq(1, max(X)), r, p), col="sienna", lwd=2)
}

### (c)
N = 1000

par(mfrow=c(1, 3))
for (p in c(0.001, 0.01, 0.5)) {
  tmean = -(1-p)/(p*log(p))
  tvar = -(1-p)/(p^2*log(p)) - tmean^2
  tstd = sqrt(tvar)
  mrange = round(seq(max(1, tmean - 3*tstd), tmean + 3*tstd, 1))
  density = function(x, p) {
    -(1-p)^x/(x*log(p))
  }
  prob = cumsum(density(mrange, p))
  X = rep(1, N)
  
  for (i in 1:N) {
    u = runif(1)
    X[i] = mrange[1] + sum(prob < u)
  }
  
  M = max(X)
  
  hist(X, freq=FALSE, main=paste("Logarithmic Series with p=", p, sep=""))
  lines(density(seq(1, M), p), col="sienna", lwd=2)
}



## 2.18

N = 10000
supp = seq(-3, 3, 0.001)

h = function(x) { exp(-x^2/2)*(sin(6*x)^2 + 3*cos(x)^2*sin(4*x)^2 + 1) }
g = dnorm
M = 5 * sqrt(2 * pi)
X = rnorm(N)
U = runif(N)
Y = M * U * g(X)

AcceptX = X[Y <= h(X)]
RejectX = X[Y > h(X)]
rate = length(AcceptX) / N
c = rate * M
AcceptY = Y[Y <= h(X)] / c
RejectY = Y[Y > h(X)] / c
f = function(x) { h(x) / c }

par(mfrow=c(1, 1))
plot(RejectX, RejectY, xlab="", ylab="", pch=20, cex=0.5, col="gray")
points(AcceptX, AcceptY, pch=20, cex=0.5)
lines(supp, g(supp) / rate, col="orange", lwd=2)
lines(supp, f(supp), col="red", lwd=2)
legend("topright", legend=c("g(x)", "f(x)"), col=c("orange", "red"), lty=c(1, 1))


## 2.19

N = 10000
supp = seq(-10, 10, 0.001)

f = dnorm
g = function(x) { exp(-abs(x))/2 }
M = sqrt(2 * exp(1) / pi)
X = rexp(N) * (2*rbinom(N, 1, 0.5) - 1)
U = runif(N)
Y = M * U * g(X)

AcceptX = X[Y <= f(X)]
RejectX = X[Y > f(X)]
rate = length(AcceptX) / N
c = rate * M
AcceptY = Y[Y <= f(X)]
RejectY = Y[Y > f(X)]

par(mfrow=c(1, 1))
plot(RejectX, RejectY, xlab="", ylab="", pch=20, cex=0.5, col="gray")
points(AcceptX, AcceptY, pch=20, cex=0.5)
lines(supp, g(supp) / rate, col="orange", lwd=2)
lines(supp, f(supp), col="red", lwd=2)
legend("topright", legend=c("g(x)", "f(x)"), col=c("orange", "red"), lty=c(1, 1))


## 2.20
### (a)

N = 500
supp = seq(-5, 5, 0.001)

f = dnorm
g = dcauchy
M = 2 / sqrt(exp(1))
X = rcauchy(N)
U = runif(N)
Y = M * U * g(X)

AcceptX = X[Y <= f(X)]
RejectX = X[Y > f(X)]
rate = length(AcceptX) / N
c = rate * M
AcceptY = Y[Y <= f(X)]
RejectY = Y[Y > f(X)]

plot(RejectX, RejectY, xlab="", ylab="", pch=20, cex=0.5, col="gray")
points(AcceptX, AcceptY, pch=20, cex=0.5)
lines(supp, g(supp) / rate, col="orange", lwd=2)
lines(supp, f(supp), col="red", lwd=2)
legend("topright", legend=c("g(x)", "f(x)"), col=c("orange", "red"), lty=c(1, 1))

### (b)

N = 2500
supp = seq(-10, 10, 0.001)

f = dnorm
g = function(x) { exp(-abs(x))/2 }
M = sqrt(2 * exp(1) / pi)
X = rexp(N) * (2*rbinom(N, 1, 0.5) - 1)
U = runif(N)
Y = M * U * g(X)

AcceptX = X[Y <= f(X)]
RejectX = X[Y > f(X)]
rate = length(AcceptX) / N
c = rate * M
AcceptY = Y[Y <= f(X)]
RejectY = Y[Y > f(X)]

par(mfrow=c(1, 1))
plot(RejectX, RejectY, xlab="", ylab="", pch=20, cex=0.5, col="gray")
points(AcceptX, AcceptY, pch=20, cex=0.5)
lines(supp, g(supp) / rate, col="orange", lwd=2)
lines(supp, f(supp), col="red", lwd=2)
legend("topright", legend=c("g(x)", "f(x)"), col=c("orange", "red"), lty=c(1, 1))