

test_models = function(){

  library(lme4)
  library(MASS)
  library(glmmTMB)
  library(fastglmm)

  object <- lmer(Reaction ~ Days + I(Days^2) + (1 | Subject), data = sleepstudy)
  varpart(object)

  # Fit a sample lme4 model
  model <- lmer(Reaction ~ Days + I(Days^2) + (Days | Subject), data = sleepstudy)
  varpart(model)

  model <- fastlmm(Reaction ~ Days + I(Days^2) + (1 | Subject), data = sleepstudy)
  varpart(model)

  gm1 <- glmer(cbind(incidence, size - incidence) ~ (1|period) + (1 | herd), data = cbpp, family = binomial)
  varpart(gm1)

  gm1 <- glmmTMB(cbind(incidence, size - incidence) ~ (1|period) + (1 | herd), data = cbpp, family = binomial)
  varpart(gm1)

  m1 <- glmmTMB(count ~ mined + (1|site),
     family=poisson, data=Salamanders)
  varpart(m1)

  m1 <- glmmTMB(count ~ mined + (Wtemp|site),
     family=nbinom2, data=Salamanders)
  varpart(m1)

  m1 <- glmmTMB(count ~ mined + (1|site/spp),
     family=nbinom2, data=Salamanders)
  varpart(m1)

  quine.nb1 <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine)
  varpart(quine.nb1)

  m1 <- glm(count ~ mined + spp,
     family=negative.binomial(1.6), data=Salamanders)
  varpart(m1)

  m1 <- lm(log(count+1) ~ mined + spp, data=Salamanders)
  varpart(m1)
}

test_getLambda = function(){

  # test approximation of lambda form count models

  # q()
  # R
  library(fastglmm)
  library(lme4)
  library(MASS)
  data(PsychAD)

  # fastglmm.nb
  form <- PTPRG ~ offset(log(libSize)) + (1|SubID)
  fit <- fastglmm.nb(form, PsychAD)
  fastglmm:::getLambda(fit)
  fastglmm:::getLambda(fit, method = "mean")
  varpart(fit)
  # varpart(fit, lambda.method = "mean")

   # fastglmm.nb
  form <- PTPRG ~ offset(log(libSize)) + Dx + (1|SubID)
  fit <- fastglmm.nb(form, PsychAD)
  fit.null <- fastglmm.nb(form, PsychAD)
  varpart(fit)
  # varpart(fit, lambda.method = "mean")


  # glm
  form <- PTPRG ~ offset(log(libSize)) + Dx
  fit <- glm(form, PsychAD, family=negative.binomial(10))
  fastglmm:::getLambda(fit)
  fastglmm:::getLambda(fit, method = "mean")
  varpart(fit)
  # varpart(fit, lambda.method = "mean")

  # glm.nb
  form <- PTPRG ~ offset(log(libSize)) + Dx
  fit <- glm.nb(form, PsychAD)
  fastglmm:::getLambda(fit)
  fastglmm:::getLambda(fit, method = "mean")
  varpart(fit)
  # varpart(fit, lambda.method = "mean")

  # library(glmmTMB)

  # form <- PTPRG ~ offset(log(libSize)) + Dx + (1|SubID)
  # fit = glmmTMB(form, data = PsychAD, family=nbinom2)
  #  varpart(fit)

}

test_varpart = function(){

  library(fastglmm)
  library(glmmTMB)
  library(MASS)
  library(mvtnorm)
  library(insight)
  library(Matrix)
  library(performance)
  library(RUnit)
  library(lme4)

  # fastglmm::varpart() gives the model R2 and 1 - [Residual frac]
  # This is equivalent to the performance::r2_nakagawa() with the trigamma method
  # Any difference should be due to parameter estimates from 
  # fastglmm() vs glmmTMB()
  # As of Jul 21, 2026 varpart uses QL dispersion

  # insight::get_variance(fit.tmb)

  set.seed(1)
  n = 1e5
  sigSq_g = 1
  beta = c(1,1)
  rho = .7
  Sigma = matrix(c(1,rho,rho,1), 2)
  X = scale(rmvnorm(n, c(0,0), Sigma))
  z = factor(sample(seq(100), n, replace=TRUE))
  df = data.frame(X, z)
  Z = sparse.model.matrix(~z, df)
  alpha = rnorm(ncol(Z), 0, sd=sqrt(sigSq_g))
  eta = as.matrix(X %*% beta + Z %*% alpha) - 1
  df$y.poisson = rpois(n, exp(eta))
  df$z.perm = factor(sample(seq(2), n, replace=TRUE))


  # weighted lm
  ##############
  w = sqrt(seq(1, nrow(df)))
  # w = w / mean(w)
  fit1 = fastlmm(log(y.poisson+1) ~ X1 + X2 + (1|z.perm), df, weights=w)
  fit2 = fastglmm(log(y.poisson+1) ~ X1 + X2 + (1|z.perm), df, weights=w)
  fit3 = lm(log(y.poisson+1) ~ X1 + X2, df, weights=w)
  fit4 = glm(log(y.poisson+1) ~ X1 + X2, df, weights=w, family=gaussian())
  fit5 = lmer(log(y.poisson+1) ~ X1 + X2 + (1|z.perm), df, weights=w)

  # fastglmm and fastlmm are slightly different
  # due to convergence
  coef(fit1)
  coef(fit2)

  fastglmm:::vpOther(fit1)
  fastglmm:::vpOther(fit2)

  # Compare to R2
  v1 = summary(fit3)$r.squared
  v2 = 1 - varpart(fit3)[3]
  checkEqualsNumeric(v1, v2, tol=1e-5)

  # sigma(fit1)
  # sigma(fit2)
  # sigma(fit3)
  # sigma(fit4)


  # fastglmm:::get_mean_weights(fit1)
  # fastglmm:::get_mean_weights(fit2)
  # fastglmm:::get_mean_weights(fit3)
  # fastglmm:::get_mean_weights(fit4)


  # fastglmm:::getDistrVar(fit1)
  # fastglmm:::getDistrVar(fit2)
  # fastglmm:::getDistrVar(fit3)
  # fastglmm:::getDistrVar(fit4)


  checkEqualsNumeric(sigma(fit1), sigma(fit2), tol=1e-3)
  checkEqualsNumeric(sigma(fit1), sigma(fit3), tol=1e-3)
  checkEqualsNumeric(sigma(fit1), sigma(fit4), tol=1e-3)
  checkEqualsNumeric(sigma(fit1), sigma(fit5), tol=1e-3)

  checkEqualsNumeric(varpart(fit1), varpart(fit2), tol=5e-2)
  checkEqualsNumeric(varpart(fit1)[-3], varpart(fit3), tol=5e-2)
  checkEqualsNumeric(varpart(fit1)[-3], varpart(fit4), tol=5e-2)
  checkEqualsNumeric(varpart(fit1), varpart(fit5), tol=5e-2)

  # Poisson model
  ###############
  fam = poisson()
  fit.tmb = glmmTMB(y.poisson ~ X1 + X2 + (1|z), df, family=fam)
  fit = fastglmm(y.poisson ~ X1 + X2 + (1|z), df,family=fam)

  # distributional variance should be very close
  get_variance_distribution( fit.tmb, approximation="trigamma" )
  fastglmm:::getDistrVar( fit, method="trigamma" )

  # Coefficient of determination
  res1 = r2_nakagawa(fit.tmb, approximation="trigamma")
  res2 = varpart(fit, method="trigamma")

  checkEqualsNumeric(1 - res1$R2_conditional,  
                  res2['Residuals'], 
                  tol = 1e-1)

  # Compare GLM with GLMM with zero variance component
  df$z.perm = factor(sample(seq(2), n, replace=TRUE))
  fit = fastglmm(y.poisson ~ X1 + X2 + (1|z.perm), df, family=fam)
  fit2 = glm(y.poisson ~ X1 + X2, df, family=fam)

  checkEqualsNumeric( varpart(fit)[-3], varpart(fit2), tol=1e-3)

  # NB model
  ##########
  df$y.nb = rnegbin(n, exp(eta), 5)
  fam = negative.binomial(NA)
  fit.tmb = glmmTMB(y.nb ~ X1 + X2 + (1|z), df, family=nbinom2)
  fit = fastglmm(y.nb ~ X1 + X2 + (1|z), df,family=fam)

  # distributional variance should be very close
  get_variance_distribution( fit.tmb, approximation="trigamma" )
  fastglmm:::getDistrVar( fit, method="trigamma" )

  # Coefficient of determination
  res1 = r2_nakagawa(fit.tmb, approximation="trigamma")
  res2 = varpart(fit, method="trigamma")

  checkEqualsNumeric(1 - res1$R2_conditional,  
                  res2['Residuals'] + res2['CountNoise'], 
                  tol = 1e-1)

  # Compare GLM with GLMM with zero variance component
  fam = negative.binomial(10)
  df$z.perm = factor(sample(seq(2), n, replace=TRUE))
  fit = fastglmm(y.nb ~ X1 + X2 + (1|z.perm), df, family=fam)
  fit2 = glm(y.nb ~ X1 + X2, df, family=fam)

  checkEqualsNumeric( varpart(fit)[-3], varpart(fit2), tol=1e-3)

  # logit
  ########
  df$y = rbinom(n, 1, plogis(eta))
  fam = binomial("logit")
  fit.tmb = glmmTMB(y ~ X1 + X2 + (1|z), df, family=fam)
  fit = fastglmm(y ~ X1 + X2 + (1|z), df,family=fam)

  # distributional variance should be very close
  get_variance_distribution( fit.tmb, approximation="trigamma" )
  fastglmm:::getDistrVar( fit, method="trigamma" )

  # Coefficient of determination
  res1 = r2_nakagawa(fit.tmb, approximation="trigamma")
  res2 = varpart(fit, method="trigamma")

  checkEqualsNumeric(1 - res1$R2_conditional,  
                  res2['Residuals'], 
                  tol = 1e-1)

  # Compare GLM with GLMM with zero variance component
  df$z.perm = factor(sample(seq(2), n, replace=TRUE))
  fit = fastglmm(y ~ X1 + X2 + (1|z.perm), df, family=fam)
  fit2 = glm(y ~ X1 + X2, df, family=fam)

  checkEqualsNumeric( varpart(fit)[-3], varpart(fit2), tol=1e-3)

  # continuous beta
  #################
  df$y = plogis(eta + rnorm(n))
  fam = binomial()
  fit.tmb = glmmTMB(y ~ X1 + X2 + (1|z), df, family=fam)
  fit = fastglmm(y ~ X1 + X2 + (1|z), df,family=fam)

  # distributional variance should be very close
  get_variance_distribution( fit.tmb, approximation="trigamma" )
  fastglmm:::getDistrVar( fit, method="trigamma" )

  # Coefficient of determination
  res1 = r2_nakagawa(fit.tmb, approximation="trigamma")
  res2 = varpart(fit, method="trigamma")

  checkEqualsNumeric(1 - res1$R2_conditional,  
                  res2['Residuals'], 
                  tol = 1e-1)

  # Compare GLM with GLMM with zero variance component
  df$z.perm = factor(sample(seq(2), n, replace=TRUE))
  fit = fastglmm(y ~ X1 + X2 + (1|z.perm), df, family=fam)
  fit2 = glm(y ~ X1 + X2, df, family=fam)

  checkEqualsNumeric( varpart(fit)[-3], varpart(fit2), tol=1e-3)

  # Compare GLM with GLMM with zero variance component
  df$z.perm = factor(sample(seq(2), n, replace=TRUE))
  fam = gaussian()
  fit = fastglmm(y ~ X1 + X2 + (1|z.perm), df, family=fam)
  fit2 = glm(y ~ X1 + X2, df, family=fam)

  checkEqualsNumeric( varpart(fit)[-3], varpart(fit2), tol=1e-3)

  # Compare lm(), glm() and r.squared
  w = seq(nrow(df))^2
  # w = w / mean(w)
  fit = lm(y ~ X1 + X2, data=df, weights=w)
  v1 = summary(fit)$r.squared
  v2 = 1 - varpart(fit)[3]

  fit2 = glm(y ~ X1 + X2, data=df, weights=w)
  v3 = 1 - varpart(fit2)[3]

  checkEqualsNumeric(v1, v2, tol = 1e-3)
  checkEqualsNumeric(v1, v3, tol = 1e-3)
}








