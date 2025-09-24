


# devtools::reload("/Users/gabrielhoffman/workspace/repos/fastglmm")


# q()
# R

# library(fastglmm)
# library(Matrix)
# library(mvtnorm)
# library(MASS)
# library(lme4)


# n = 1e6 
# set.seed(1)
# rho = .999
# beta = c(1,2)
# sigSq_g = 1e-10
# set.seed(1)
# Sigma = matrix(c(1,rho,rho,1), 2)
# X = scale(rmvnorm(n, c(0,0), Sigma))
# z = factor(sample(seq(2), n, replace=TRUE))
# df = data.frame(X, z)
# Z = sparse.model.matrix(~z, df)
# alpha = rnorm(ncol(Z), 0, sd=sqrt(sigSq_g))
# eta = 10 + as.matrix(X %*% beta + Z %*% alpha)

# df$y.logistic = rbinom(n, 1, prob=plogis(eta))
# df$y.poisson = rpois(n, exp(eta))
# df$y.nb = MASS::rnegbin(n, exp(eta), theta = 100000)

# fam = quasipoisson()

# form = y.nb ~ X1 + X2
# fit.glm = glm(form, df, family=fam)

# dispersion(fit.glm)

# form = y.nb ~ X1 + X2 + (1|z)
# fit.fast = fastglmm(form, df, family=fam)


# dispersion(fit.fast)



# form = y.nb ~ X1 + X2 
# fit.pql = glmmPQL(form, random = ~ 1|z, df, family=fam)

# coef(summary(fit.fast))
# coef(summary(fit.pql))


# dispersion(fit.pql)
# dispersion(fit.fast)

# vcov(fit.glm)
# vcov(fit.pql)
# vcov(fit.fast)





# dispersion(fit.fast)
# summary(fit.fast)$dispersion


# vcov(fit.glm)
# vcov(fit.fast)

# summary(fit.glm)
# summary(fit.fast)