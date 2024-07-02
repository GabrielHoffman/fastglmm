

Z = Rfast::matrnorm(4000, 3000)
system.time(dcmp2 <- svd(Z))

library(MASS)
library(lme4)
library(glmmTMB)
library(fastglmm)
set.seed(1)

n = 300000
ndonors = 1500

info = data.frame(x = rnorm(n))
info$Indiv = factor(sample(seq(ndonors), n, replace=TRUE))
info$Indiv = droplevels(info$Indiv)
beta = .0
eta = 0 + info$x * beta + model.matrix(~ 0 + Indiv, info) %*% rnorm(nlevels(info$Indiv))
info$y = rnegbin(n, mu=exp(eta), theta = 10)
# hist(info$y)

# Poisson
fit1 <- glmer(y ~ x + (1|Indiv), data = info, family="poisson")
coef(summary(fit1))

# NB
fit2 <- glmer.nb(y ~ x + (1|Indiv), data = info)
coef(summary(fit2))

fit2a <- glmmTMB(y ~ x + (1|Indiv), data=info, family=nbinom2)
coef(summary(fit2a))

# quasi quadratic
fam <- quasi(link = "log", variance="mu^2")
fit3 = glmmPQL(y ~ x, random = ~ 1|Indiv, data = info, family=fam, verbose=FALSE)
coef(summary(fit3))

# quasi linear
# estimates are same as Poisson by definition
fam <- quasi(link = "log", variance="mu")
fit4 = glmmPQL(y ~ x, random = ~ 1|Indiv, data = info, family=fam, verbose=FALSE)
coef(summary(fit4))






n_reps = 1
system.time(replicate(n_reps, 
	glmer(y ~ x + (1|Indiv), data = info, family="poisson")))

system.time(replicate(n_reps, 
	glmer.nb(y ~ x + (1|Indiv), data = info)))

system.time(replicate(n_reps, 
	glmmTMB(y ~ x + (1|Indiv), data=info, family=nbinom2)))

system.time(replicate(n_reps, 
	glmmPQL(y ~ x, random = ~ 1|Indiv, data = info, family=fam, verbose=FALSE)))

# glmer.nb() fits glmer() twice
# glmmPQL() evalutes a sequence of reweighted nlme() calls 



# Poisson
fit = glm(y ~ x, family="poisson")
coef(summary(fit))

# NB
fit = glm.nb(y ~ x)
coef(summary(fit))

# quasi
fit = glm(y ~ x, family=quasi(link = "log", variance="mu^2"))
coef(summary(fit))


# lrgpr
######
library(lrgpr)
library(Matrix)

# Compute SVD directly from info$Indiv, using count order
# dsgn = sparse.model.matrix(~ 0 + Indiv, info)
Z = t(fac2sparse(info$Indiv))
ZtZ = crossprod(Z)
dcmp = eigen(ZtZ, symmetric=TRUE)
# Z = UDV^T
# Z V^-T D^-1 = U
# ZtZ = V D U^T U D V^T
# dcmp2 = svd(Z)

dcmp$vectors = tcrossprod(Z, solve(dcmp$vectors)) %*% diag(1/sqrt(dcmp$values))
# dcmp$vectors = Matrix(dcmp$vectors, sparse = TRUE)
dcmp$vectors = as.matrix(dcmp$vectors)

U[1:2,]
dcmp2$u[1:2,]

A = crossprod(U, dcmp2$u)
table(as.numeric(zapsmall(A)))


fit = lrgpr( info$y ~ info$x, dcmp)

X = model.matrix(~ x, info)
# source("/Users/gabrielhoffman/workspace/repos/fastlmm/R/fastlmm.R")
fit = fastlmm_R( info$y, X, U = dcmp$vectors, s = sqrt(dcmp$values))
fit$beta


X = model.matrix(~x, info)
fit = fastlmm_R( info$y, X, U = dcmp$vectors, s = sqrt(dcmp$values))

fastglmm:::fastlmm_c(info$y, X, U = dcmp$vectors, s = sqrt(dcmp$values))


# give same answer
fit1 <- lmer(y ~ x + (1|Indiv), data = info, REML=FALSE)
coef(summary(fit1))



system.time(replicate(10, 
	lrgpr( info$y ~ info$x, dcmp)))

system.time(replicate(100, 
	lrgpr( info$y ~ info$x, dcmp, delta=0.5299353)))


system.time(replicate(500, 
	fastlmm( info$y, X, U = dcmp$vectors, s = sqrt(dcmp$values))))


system.time(replicate(500, 
	fastlmm(info$y, X, U = dcmp$vectors, s = sqrt(dcmp$values), delta = 0.008580024)))


system.time(replicate(100, 
	lmer(y ~ x + (1|Indiv), data = info, REML=FALSE)))

fit = lmer(y ~ x + (1|Indiv), data = info, REML=FALSE)
system.time(replicate(100, 
	refit(fit, newweights = runif(nrow(info)))
	))



info$y = as.numeric(info$y)
library(profvis)
profvis({
	replicate(100, 
	fastlmm( info$y, X, U = dcmp$vectors, s = sqrt(dcmp$values)))
	})


library(profvis)
profvis({
	replicate(100, 
	fastlmm( info$y, X, U = dcmp$vectors, s = sqrt(dcmp$values), delta=1))
	})






library(profvis)
profvis({
	replicate(1000, 
	lmer(y ~ x + (1|Indiv), data = info, REML=FALSE))
	})



object <- glmer(y ~ x + (1|Indiv), data = info, family="poisson")

resp <- model.response(model.frame(object))
mu <- na.omit(fitted(object))
system.time(
res1 <- MASS::theta.ml(resp, mu, weights = object@resp$weights, limit=20,
				eps = .Machine$double.eps^0.5)
)


# devtools::reload("/Users/gabrielhoffman/workspace/repos/fastglmm")
system.time(
res2 <- fastglmm:::theta_ml(resp, mu, n = length(mu),
				weights = object@resp$weights,
				limit = 20,
				eps = .Machine$double.eps^0.5)
)


res1
res2



