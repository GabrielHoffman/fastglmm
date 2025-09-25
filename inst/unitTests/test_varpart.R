

test_varpart = function(){

  library(fastglmm)
  library(glmmTMB)
  library(MASS)
  library(mvtnorm)
  library(insight)
  library(performance)
  library(RUnit)

  # fastglmm::varpart() gives the model R2 and 1 - [Residual frac]
  # This is equivalent to the performance::r2_nakagawa() with the trigamma method
  # Any difference should be due to parameter estimates from 
  # fastglmm() vs glmmTMB()

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
  eta = as.matrix(X %*% beta + Z %*% alpha) - 4
  df$y.poisson = rpois(n, exp(eta))

  # Poisson model
  fam = poisson()
  fit.tmb = glmmTMB(y.poisson ~ X1 + X2 + (1|z), df, family=fam)
  fit.null = glmmTMB(y.poisson ~ (1|z), df, family=fam)
  fit = fastglmm(y.poisson ~ X1 + X2 + (1|z), df,family=fam)

  # distributional variance should be very close
  get_variance_distribution( fit.tmb, null_model = fit.null, approximation="trigamma" )
  fastglmm:::getDistrVar( fit, method="trigamma" )

  # Coefficient of determination
  res1 = r2_nakagawa(fit.tmb, null_model = fit.null, approximation="trigamma")
  res2 = varpart(fit)

  checkEqualsNumeric(1 - res1$R2_conditional,  
                  res2['Residuals'], 
                  tol = 1e-1)


  # NB model
  df$y.nb = rnegbin(n, exp(eta), 5)
  fam = negative.binomial(NA)
  fit.tmb = glmmTMB(y.nb ~ X1 + X2 + (1|z), df, family=nbinom2)
  fit = fastglmm(y.nb ~ X1 + X2 + (1|z), df,family=fam)

  # distributional variance should be very close
  get_variance_distribution( fit.tmb, approximation="trigamma" )
  fastglmm:::getDistrVar( fit, method="trigamma" )

  # Coefficient of determination
  res1 = r2_nakagawa(fit.tmb, approximation="trigamma")
  res2 = varpart(fit)

  checkEqualsNumeric(1 - res1$R2_conditional,  
                  res2['Residuals'], 
                  tol = 1e-1)


  # logit
  df$y = rbinom(n, 1, plogis(eta))
  fam = binomial("logit")
  fit.tmb = glmmTMB(y ~ X1 + X2 + (1|z), df, family=fam)
  fit = fastglmm(y ~ X1 + X2 + (1|z), df,family=fam)

  # distributional variance should be very close
  get_variance_distribution( fit.tmb, approximation="trigamma" )
  fastglmm:::getDistrVar( fit, method="trigamma" )

  # Coefficient of determination
  res1 = r2_nakagawa(fit.tmb, approximation="trigamma")
  res2 = varpart(fit)

  checkEqualsNumeric(1 - res1$R2_conditional,  
                  res2['Residuals'], 
                  tol = 1e-1)


  # continuous beta
  df$y = plogis(eta + rnorm(n))
  fam = binomial()
  fit.tmb = glmmTMB(y ~ X1 + X2 + (1|z), df, family=fam)
  fit = fastglmm(y ~ X1 + X2 + (1|z), df,family=fam)

  # distributional variance should be very close
  get_variance_distribution( fit.tmb, approximation="trigamma" )
  fastglmm:::getDistrVar( fit, method="trigamma" )

  # Coefficient of determination
  res1 = r2_nakagawa(fit.tmb, approximation="trigamma")
  res2 = varpart(fit)

  checkEqualsNumeric(1 - res1$R2_conditional,  
                  res2['Residuals'], 
                  tol = 1e-1)

}
