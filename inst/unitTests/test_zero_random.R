


test_zero_random = function(){
  library(fastglmm)
  library(MASS)
  library(RUnit)

  set.seed(3)

  # Simulate large dataset where random effect is zero, sigSq_g is set to small value and then compared to standard GLM
  n = 10000
  data = data.frame(Age = runif(n, 0, 100))
  data$ID = factor(sample(seq(1000), n, replace=TRUE))

  design = model.matrix( ~ Age, data)
  eta = design %*% c(2, .01)
  data$y.count = rpois(n, exp(eta))

  form1 = y.count ~ scale(Age) + (1|ID)
  form2 = y.count ~ scale(Age) 

  # Poisson
  fam = poisson() 
  fit1 = fastglmm( form1, data, family=fam, delta=1e8)
  fit2 = glm( form2, data, family=fam )

  checkEqualsNumeric(coef(fit1), coef(fit2), tol=1e-5)
  checkEqualsNumeric(vcov(fit1), vcov(fit2), tol=1e-4)
  checkEqualsNumeric(dispersion(fit1), 
    summary(fit2)$dispersion, tol=1e-5)
  checkEqualsNumeric(residuals(fit1, "pearson"), 
                      residuals(fit2, "pearson"), tol=1e-5)
  checkEqualsNumeric(residuals(fit1, "deviance"), 
                      residuals(fit2, "deviance"), tol=1e-5)


  # Quasipoisson
  data$y.count = rnegbin(n, exp(eta), 10)
  fam = quasipoisson() 
  fit1 = fastglmm( form1, data, family=fam, delta=1e8)
  fit2 = glm( form2, data, family=fam )

  dispersion(fit1)

  checkEqualsNumeric(coef(fit1), coef(fit2), tol=1e-5)
  checkEqualsNumeric(vcov(fit1), vcov(fit2), tol=1e-4)
  checkEqualsNumeric(dispersion(fit1), 
    summary(fit2)$dispersion, tol=1e-5)
  checkEqualsNumeric(residuals(fit1, "pearson"), 
                      residuals(fit2, "pearson"), tol=1e-5)
  checkEqualsNumeric(residuals(fit1, "deviance"), 
                      residuals(fit2, "deviance"), tol=1e-5)


  # Binomial
  eta = design %*% c(1, .01) + rnorm(nrow(design))
  data$y01 = rbinom(n, 1, plogis(eta))

  form1 = y01 ~ scale(Age) + (1|ID)
  form2 = y01 ~ scale(Age) 

  fam = binomial() 
  fit1 = fastglmm( form1, data, family=fam, delta=1e8)
  fit2 = glm( form2, data, family=fam )

  checkEqualsNumeric(coef(fit1), coef(fit2), tol=1e-5)
  checkEqualsNumeric(vcov(fit1), vcov(fit2), tol=1e-3)
  checkEqualsNumeric(dispersion(fit1), 
    summary(fit2)$dispersion, tol=1e-5)
  checkEqualsNumeric(residuals(fit1, "pearson"), 
                      residuals(fit2, "pearson"), tol=1e-5)
  checkEqualsNumeric(residuals(fit1, "deviance"), 
                      residuals(fit2, "deviance"), tol=1e-5)



  # Quasibinomial
  eta = design %*% c(1, .01) + rnorm(nrow(design))
  data$beta = plogis(eta)

  form1 = beta ~ scale(Age) + (1|ID)
  form2 = beta ~ scale(Age) 

  fam = quasibinomial() 
  fit1 = fastglmm( form1, data, family=fam, delta=1e8)
  fit2 = glm( form2, data, family=fam )

  checkEqualsNumeric(coef(fit1), coef(fit2), tol=1e-5)
  checkEqualsNumeric(vcov(fit1), vcov(fit2), tol=1e-4)
  checkEqualsNumeric(dispersion(fit1), 
    summary(fit2)$dispersion, tol=1e-5)
  checkEqualsNumeric(residuals(fit1, "pearson"), 
                      residuals(fit2, "pearson"), tol=1e-5)
  checkEqualsNumeric(residuals(fit1, "deviance"), 
                      residuals(fit2, "deviance"), tol=1e-5)
}
