# 


test_quasidispersion = function(){
  
  library(MASS)
  library(fastglmm)
  library(RUnit)

  set.seed(101)

  dd <- expand.grid(
    f1 = factor(1:4),
    f2 = LETTERS[1:2], 
    g = factor(1:12), 
    rep = 1:15,
    KEEP.OUT.ATTRS=FALSE)

  mu <- 5*(-4 + with(dd, as.integer(f1) + 4*as.numeric(f2)))
  mu = mu + model.matrix(~0+g, dd) %*% rnorm(12, sd=1e-12)
  dd$y <- rnbinom(nrow(dd), mu = mu, size = 0.5)

  Y = matrix(as.numeric(dd$y), nrow=1)
  design = model.matrix(y ~ f1*f2 , dd)
  Z = preprocess_indicator(dd$g)

  # NB estimate theta
  ####################
  fit1 = glm.nb(y ~ f1*f2, data = dd)
  fit3 = fastglmm(y ~ f1*f2 + (1|g), data =dd, family=negative.binomial(NA))

  a = summary(fit1)$dispersion
  b = fit3$dispersion
  checkEqualsNumeric(a,b)

  a = sqrt(diag(vcov(fit1)))
  b = sqrt(diag(vcov(fit3)))
  checkEqualsNumeric(a,b, tol=1e-2)

  # NB fixed theta
  ################
  fam = negative.binomial(52)
  fit1 = glm(y ~ f1*f2, data = dd, family = fam)
  fit3 = fastglmm(y ~ f1*f2 + (1|g), data =dd, family=fam)

  a = summary(fit1)$dispersion
  b = fit3$dispersion
  checkEqualsNumeric(a,b, tol=1e-3)

  a = sqrt(diag(vcov(fit1)))
  b = sqrt(diag(vcov(fit3)))
  checkEqualsNumeric(a,b, tol=1e-2)

  # Poisson
  #########

  fam = poisson()
  fit1 = glm(y ~ f1*f2, data = dd, family = fam)
  fit3 = fastglmm(y ~ f1*f2 + (1|g), data =dd, family=fam)

  a = summary(fit1)$dispersion
  b = fit3$dispersion
  checkEqualsNumeric(a,b, tol=1e-3)

  a = sqrt(diag(vcov(fit1)))
  b = sqrt(diag(vcov(fit3)))
  checkEqualsNumeric(a,b, tol=1e-2)

 # Quasi-Poisson 
  ##############

  fam = quasipoisson()
  fit1 = glm(y ~ f1*f2, data = dd, family = fam)
  fit3 = fastglmm(y ~ f1*f2 + (1|g), data =dd, family=fam)

  a = summary(fit1)$dispersion
  b = fit3$dispersion
  checkEqualsNumeric(a,b, tol=1e-3)

  a = sqrt(diag(vcov(fit1)))
  b = sqrt(diag(vcov(fit3)))
  checkEqualsNumeric(a,b, tol=1e-2)

}







  