
test_ddf = function(){

  library(fastglmm)
  library(lmerTest)
  library(RUnit)

  set.seed(1)
  dd <- expand.grid(
    g = factor(1:10), 
    rep = 1:10,
    KEEP.OUT.ATTRS=FALSE)

  dd$x <- rnorm(nrow(dd))
  mu <- 5*dd$x
  mu <- mu + model.matrix(~0+g, dd) %*% rnorm(nlevels(dd$g), sd=1)

  dd$y <- mu + rnorm(nrow(mu))

  # fastlmm
  fit <- fastlmm(y ~ x + (1 | g), dd)
  ddf_fastlmm = fastglmm:::ddf(fit)

  # Compare to lmerTest
  #####################

  fit2 <- lmer(y ~ x + (1 | g), dd, REML=FALSE)
  ddf_lmertest <- coef(summary(fit2))[,3]

  checkEqualsNumeric(ddf_fastlmm, ddf_lmertest, tol=1e-3)
  checkEqualsNumeric(coef(summary(fit)), coef(summary(fit2)), tol=1e-3)

  # Compare interfaces
  #####################
  a = coef(summary(fit))
  b = anova(fit)

  checkEqualsNumeric(a[,"ddf"], b["df2"], tol=1e-3)
  checkEqualsNumeric(log(a[,'Pr(>|t|)']), log(b[,4]))

  d = linearHypothesis(fit, c(1,0), test="F")
  checkEqualsNumeric(d$`Pr(>F)`[2], b[1,4])

}


test_ddf_glmmTMB = function(){

  library(fastglmm)
  library(glmmTMB)
  library(RUnit)

  set.seed(1)
  dd <- expand.grid(
    g = factor(1:10), 
    rep = 1:10,
    KEEP.OUT.ATTRS=FALSE)

  dd$x <- rnorm(nrow(dd))
  mu <- dd$x
  mu <- 4 + mu + model.matrix(~0+g, dd) %*% rnorm(nlevels(dd$g), sd=1)

  dd$y <- rpois(nrow(dd), exp(mu))

  # Compare to glmmTMB: gaussian
  ################################

  family = gaussian()
  fit <- fastglmm(y ~ x + (1 | g), dd, family=family)
  ddf_fastglmm = fastglmm:::ddf(fit)

  fit2 <- glmmTMB(y ~ x + (1 | g), dd, family=family)
  ddf_tmb <- coef(summary(fit2, ddf="satterthwaite"))$cond[,"ddf"]

  checkEqualsNumeric(ddf_fastglmm, ddf_tmb, tol=1e-2)

  res1 = coef(summary(fit))
  res2 = coef(summary(fit2, ddf="satterthwaite"))$cond
  checkEqualsNumeric(res1, res2[,colnames(res1)], 
    tol=1e-3)

  # Compare interfaces
  #####################

  family = poisson()
  fit <- fastglmm(y ~ x + (1 | g), dd, family=family)

  a = coef(summary(fit))
  b = anova(fit)

  checkEqualsNumeric(a[,"ddf"], b["df2"], tol=1e-3)
  checkEqualsNumeric(log(a[,'Pr(>|t|)']), log(b[,4]))

  d = linearHypothesis(fit, c(1,0), test="F")
  checkEqualsNumeric(d$`Pr(>F)`[2], b[1,4])

}



# test_ddf_glmmTMB = function(){

#   library(fastglmm)
#   library(glmmTMB)
#   library(RUnit)

#   set.seed(1)
#   dd <- expand.grid(
#     g = factor(1:102), 
#     rep = 1:10,
#     KEEP.OUT.ATTRS=FALSE)

#   dd$x <- rnorm(nrow(dd))
#   mu <- dd$x
#   mu <- 4 + mu + model.matrix(~0+g, dd) %*% rnorm(nlevels(dd$g), sd=1)

#   dd$y <- rpois(nrow(dd), exp(mu))


#   family = poisson()
#   fit <- fastglmm(y ~ x + (1 | g), dd, family=family)
#   ddf_fastglmm = fastglmm:::ddf(fit)

#   fit2 <- glmmTMB(y ~ x + (1 | g), dd, family=family)
#   ddf_tmb <- coef(summary(fit2, ddf="satterthwaite"))$cond[,"ddf"]

#   ddf_fastglmm
#   ddf_tmb


#   coef(summary(fit))
#   coef(summary(fit2))$cond


#   checkEqualsNumeric(ddf_fastglmm, ddf_tmb, tol=1e-2)

#   fit <- fastglmm(y ~ x + (1 | g), dd, family=gaussian())
#   fastglmm:::ddf(fit)

#   fit <- fastglmm(y ~ x + (1 | g), dd, family=poisson())
#   fastglmm:::ddf(fit)


#   fit2 <- glmmTMB(y ~ x + (1 | g), dd, family=gaussian())
#   coef(summary(fit2, ddf="satterthwaite"))$cond

#   fit2 <- glmmTMB(y ~ x + (1 | g), dd, family=poisson())
#   coef(summary(fit2, ddf="satterthwaite"))$cond
# }

