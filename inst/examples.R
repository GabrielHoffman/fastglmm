

library(lrgpr)
source("fastlmm/R/fastlmm.R")

set.seed(1)
n <- 5000
y <- rnorm(n)
X_geno <- matrix(rnorm(n*1000), nrow=n) 
age <- rpois(n, 50)
sex <- as.factor(sample(1:2, n, replace=TRUE))
X = model.matrix(~ sex + age)
X1 = model.matrix(~ sex -1)


# decomp <- eigen(tcrossprod(X_geno))
# fit <- lrgpr( y ~ sex + age, decomp, diagnostic=TRUE)
# fit2 = fastlmm3(y, X, decomp$vectors, decomp$values)

dcmp = svd(X_geno)
fit <- lrgpr( y ~ sex + age, dcmp, rank = 10)
fit2 = fastlmm3(y, X, dcmp$u, dcmp$d^2, rank = 10)
fit4 = fastlmm4(y, X, dcmp$u, dcmp$d^2, rank = 10)

coef(fit)
fit4$beta

logLik(fit)
fit4$ML

fit$sigSq_a
fit4$vg 

fit$sigSq_e
fit4$ve

fit$delta
fit4$delta


fit$df
fit4$df

library(microbenchmark)

n.rank = min(5000, ncol(dcmp$u))

decmp2 = list(u = dcmp$u[,seq(n.rank)], s.sq = dcmp$d[seq(n.rank)]^2)

df = microbenchmark(
	lm = lm(y ~ sex + age),
	lm.fit = lm.fit(X,y),
	lrgpr = lrgpr( y ~ sex + age, dcmp, rank = n.rank), 
	fastlmm4 = fastlmm4(y, X, dcmp$u, dcmp$d^2, rank = n.rank),  
	fastlmm4.grid10 = fastlmm4(y, X, dcmp$u, dcmp$d^2, rank = n.rank, n.grid=10),
	lrgpr.delta = lrgpr( y ~ sex + age, dcmp, rank = n.rank, delta=1),
	fastlmm4.delta = fastlmm4(y, X, dcmp$u, dcmp$d^2, rank = n.rank, delta=1), 
	times=100)



df = microbenchmark(
	fastlmm4 = fastlmm4(y, X, dcmp$u, dcmp$d^2, rank = n.rank), 
	fastlmm1 = fastlmm4(y, X1, dcmp$u, dcmp$d^2, rank = n.rank), 
	fastlmm4.decmp = fastlmm4(y, X, decmp2$u, decmp2$s.sq, rank = n.rank),
	times=10)

Rprof()
res = replicate( 1000, fastlmm4(y, X, decmp2$u, decmp2$s.sq, rank = n.rank))
res = summaryRprof()














