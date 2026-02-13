

test_refitModel = function(){
  

  # q()
  # R

  library(fastglmm)
  library(lme4)
  library(MASS)
  library(RUnit)

  # Check that refit models are identical to original

  # fastlmm
  #########

  # MLE of delta
  fit <- fastlmm(Reaction ~ Days + (1 | Subject), sleepstudy)
  fit2 <- refitModel(fit)

  sapply( names(fit), function(x){
    cat(x, "\n")
    checkEquals(fit[[x]], fit2[[x]])
    })

  # fixed delta
  fit1 <- fastlmm(Reaction ~ Days + (1 | Subject), sleepstudy, delta=1)
  fit2 <- refitModel(fit1, delta=1)

  sapply( names(fit)[-1], function(x){
    cat(x, "\n")
    checkEquals(fit1[[x]], fit2[[x]])
    })

  # interceptOnly
  fit1 <- fastlmm(Reaction ~ (1 | Subject), sleepstudy)
  fit2 <- refitModel(fit1, interceptOnly = TRUE)

  sapply( names(fit)[-1], function(x){
    cat(x, "\n")
    checkEquals(fit1[[x]], fit2[[x]])
    })

  # fastglmm
  ##########

  # MLE of delta
  fit1 = fastglmm(y ~ trt + I(week > 2) + (1 | ID),
         family = binomial(), data = bacteria)
  fit2 <- refitModel(fit1)

  sapply( names(fit1), function(x){
    cat(x, "\n")
    checkEquals(fit1[[x]], fit2[[x]])
    })

  # fixed delta
  fit1 <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
         family = binomial(), data = bacteria, delta=1)
  fit2 <- refitModel(fit1, delta=1)

  sapply( names(fit)[-1], function(x){
    cat(x, "\n")
    checkEquals(fit1[[x]], fit2[[x]])
    })

  # interceptOnly
  fit1 <- fastglmm(y ~ (1 | ID),
         family = binomial(), data = bacteria)
  fit2 <- refitModel(fit1, interceptOnly=TRUE)

  sapply( names(fit)[-1], function(x){
    cat(x, "\n")
    checkEquals(fit1[[x]], fit2[[x]])
    })



  # Negative binomial
  set.seed(1)
  n = 1e5
  sigSq_g = 1
  beta = c(1,1)
  rho = .7
  X = scale(matrix(rnorm(2*n), nrow=n))
  z = factor(sample(seq(100), n, replace=TRUE))
  df = data.frame(X, z)
  Z = sparse.model.matrix(~z, df)
  alpha = rnorm(ncol(Z), 0, sd=sqrt(sigSq_g))
  eta = as.matrix(X %*% beta + Z %*% alpha) - 4
  df$y.nb = rnegbin(n, exp(eta), 5)

  # full model
  fam = negative.binomial(NA)
  fit = fastglmm(y.nb ~ X1 + X2 + (1|z), df, family=fam, doCoxReid=FALSE)

  # null model
  fit_null = fastglmm(y.nb ~ (1|z), df, family=fam)
  fit_null2 = refitModel(fit, interceptOnly=TRUE)

  family(fit_null)$dispersion
  family(fit_null2)$dispersion


  fit_null2$data = fit_null$data = NULL
  # summary(fit_null)
  # summary(fit_null2)

  sapply( names(fit_null)[-1], function(x){
    cat(x, "\n")
    checkEquals(fit_null[[x]], fit_null2[[x]])
    })

  # Include offset term
  #####################

  eta = as.matrix(X %*% beta + Z %*% alpha) 
  df$y.nb = rnegbin(n, exp(eta), 10)

  df$off = rpois(nrow(df), df$y.nb)
  form = y.nb ~ X1 + X2 + (1|z) + offset(log(off+1))
  fit = fastglmm(form, df, family=fam)
  
  # null model
  form = y.nb ~ (1|z) + offset(log(off+1))
  fit_null = fastglmm(form, df,family=fam)
  fit_null2 = refitModel(fit, interceptOnly=TRUE)

  fit_null2$data = fit_null$data = NULL
  # summary(fit_null)
  # summary(fit_null2)

  ids = names(fit_null)
  ids = ids[!ids %in% c("formula")]
  sapply( ids, function(x){
    cat(x, "\n")
    checkEquals(fit_null[[x]], fit_null2[[x]])
    })

  fit_null$formula
  fit_null2$formula

  # Compare to GLM with offset
  ############################

  # set sigSq_g to zero
  form = y.nb ~ (1|z) + offset(log(off+1))
  fam = negative.binomial(10)
  fit = fastglmm(form, df, family=fam, delta=1e8)

  # GLM 
  form = y.nb ~ offset(log(off+1))
  fit.glm = glm(form, df, family=fam)
  
  a = fastglmm:::getLambda(fit, delta = 1e8, fixedNBtheta=TRUE)
  b = fastglmm:::getLambda(fit.glm)

  checkEqualsNumeric(a,b, tol=1e-6)
}