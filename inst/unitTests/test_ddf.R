
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
  ddf_fastlmm = ddf(fit)

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

  checkEqualsNumeric(a[,"df"], b["df2"], tol=1e-3)
  checkEqualsNumeric(log(a[,'Pr(>|t|)']), log(b[,4]))

  d = linearHypothesis(fit, c(1,0))
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
  ddf_fastglmm = ddf(fit)

  fit2 <- glmmTMB(y ~ x + (1 | g), dd, family=family)
  ddf_tmb <- coef(summary(fit2, ddf="satterthwaite"))$cond[,"ddf"]

  checkEqualsNumeric(ddf_fastglmm, ddf_tmb, tol=1e-2)

  res1 = coef(summary(fit))
  res2 = coef(summary(fit2, ddf="satterthwaite"))$cond
  checkEqualsNumeric(res1, res2[,c(1,2,4,3,5)], 
    tol=1e-3)

  # Poisson
  family = poisson()
  fit <- fastglmm(y ~ x + (1 | g), dd, family=family)
  ddf_fastglmm = ddf(fit)

  fit2 <- glmmTMB(y ~ x + (1 | g), dd, family=family)
  ddf_tmb <- coef(summary(fit2, ddf="satterthwaite"))$cond[,"ddf"]

  ddf_fastglmm
  ddf_tmb

  dd$z = fit$y
  it = lmerTest::lmer(z ~ x + (1 | g), dd, weights = fit$weights)
  coef(summary(it))


  # checkEqualsNumeric(ddf_fastglmm, ddf_tmb, tol=1e-2)

  # res1 = coef(summary(fit))
  # res2 = coef(summary(fit2, ddf="satterthwaite"))$cond
  # checkEqualsNumeric(res1, res2[,colnames(res1)], 
  #   tol=1e-3)

  # Compare interfaces
  #####################

  family = poisson()
  fit <- fastglmm(y ~ x + (1 | g), dd, family=family)

  a = coef(summary(fit))
  b = anova(fit)

  checkEqualsNumeric(a[,"df"], b["df2"], tol=1e-3)
  checkEqualsNumeric(log(a[,'Pr(>|t|)']), log(b[,4]))

  d = linearHypothesis(fit, c(1,0))
  checkEqualsNumeric(d$`Pr(>F)`[2], b[1,4])

}


test_hypothesisTesting = function(){

  library(fastglmm)
  library(lmerTest)
  library(RUnit)

  sleepstudy$Days = factor(sleepstudy$Days)

  fit1 <- lmer( Reaction ~ Days + (1 | Subject), sleepstudy, REML=FALSE)
  fit2 <- fastlmm( Reaction ~ Days + (1 | Subject), sleepstudy)

  # summary
  checkEqualsNumeric(
    coef(summary(fit1)), 
    coef(summary(fit2)),
    tol = 1e-4)

  # anova
  a = anova(fit1)
  b = anova(fit2)

  aa = as.matrix(a) 
  bb = as.matrix(b[2,])

  checkEqualsNumeric(
    aa[-c(1:2)],
    bb, 
    tol = 1e-5)

  # linearHypothesis
  res1 = linearHypothesis(fit2, "(Intercept)")

  # single variable
  checkEqualsNumeric(
    as.matrix(res1[2,]),
    as.matrix(b[1,])[c(2,1,3,4)]
  )

  # joint test
  res2 = linearHypothesis(fit2, paste0("Days", seq(1, 9)))

  checkEqualsNumeric(
    as.matrix(res2[2,]),
    as.matrix(b[2,])[c(2,1,3,4)])

}

test_ddf_compare = function(){

  # https://chatgpt.com/share/698f66b1-88b0-800b-9001-e28f444f6015
  est_hessian <- function(fit){

    delta <- fit$delta
    sigSq_g <- fit$sigSq_g
    s <- fit$s

    n <- length(fit$y)
    r <- length(s)
    H <- matrix(0, 2,2)

    a <- sum(1/(s + delta)) + (n-r)/delta
    b <- sum(1/(s + delta)^2) + (n-r)/delta^2

    H[1,1] <- (n - 2*delta*a + delta^2*b) / (2*sigSq_g^2)
    H[1,2] <- H[2,1] <- (a - delta*b) / (2*sigSq_g^2)
    H[2,2] <- b / (2*sigSq_g^2)
    H
  }

  est_gradient <- function(fit, L){

    # slow versions
    # W <- with(fit, solve(tcrossprod(Z) + diag(delta, nrow(Z))))
    # A <- crossprod(X, W) %*% X
    # B <- crossprod(X, W %*% W) %*% X

    X = diag(c(sqrt(fit$weights))) %*% fit$design
    Xu <- crossprod(fit$U, X)
    Gamma_XX <- crossprod(X) - crossprod(Xu)
    inv_s_delta <- 1 / (fit$s + fit$delta)
    inv_s_delta_Xu <- inv_s_delta * Xu

    A <- crossprod(Xu, inv_s_delta_Xu) + Gamma_XX / fit$delta

    inv_s_delta_Xu <- inv_s_delta^2 * Xu
    B <- crossprod(Xu, inv_s_delta_Xu) + Gamma_XX / fit$delta^2

    lapply(seq(nrow(L)), function(i){

      g <- c(0, 0)
      invAL <- solve(A, L[i,])
      C <- invAL %*% B %*% invAL
      g[1] <- crossprod(L[i,], invAL) - fit$delta*C
      g[2] <- C
      g
    })
  }


  ddf.orig <- function(fit, L = diag(1, length(coef(fit))) ){

    stopifnot(is(fit, "fastlmm"))

    H <- est_hessian(fit)
    g <- est_gradient(fit, L)

    sapply(seq(nrow(L)), function(i){
      var_Lbeta <- crossprod(L[i,], vcov(fit)) %*% L[i,]
      v_numerator <- 2 * var_Lbeta^2
      v_denom <- crossprod(g[[i]], solve(H)) %*% g[[i]]

      v_numerator / v_denom  
    })
  }

  library(MASS)
  library(fastglmm)
  library(Matrix)
  library(lme4)
  library(RUnit)
  set.seed(1)

  n = 100
  ndonors = 100

  # n = 3000
  # info = data.frame(x = rnorm(n))
  # info$Indiv = factor(seq(n))

  info = data.frame(x = rnorm(n))
  info$Indiv = factor(sample(seq(ndonors), n, replace=TRUE))
  info$Indiv = droplevels(info$Indiv)
  beta = 1
  eta = 4 + info$x * beta + model.matrix(~ 0 + Indiv, info) %*% rnorm(nlevels(info$Indiv), 0, sqrt(3)) 
  # info$y = rnegbin(n, mu=exp(eta), theta = 10)
  info$y = eta + rnorm(length(eta))

  dcmp = indicator_decomp( info$Indiv )
  indicObj = preprocess_indicator( info$Indiv )
  X = model.matrix( ~ x, info)

  formula = y ~ x + (1|Indiv)
  data = info
  eta = X %*% c(2, .0007) + model.matrix(~Indiv + 0, data) %*% rnorm(nlevels(info$Indiv))
  data$y = rpois(nrow(data), exp(eta))
  family = poisson()



  fit = fastglmm( y ~ x + (1|Indiv), data, family = poisson())
  a = ddf(fit)
  b = ddf.orig(fit)
  checkEqualsNumeric(a,b, tol=1e-3)


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

