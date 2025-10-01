
#' Fit generalized linear mixed model via PQL
#' 
#' Fit generalized linear mixed model (GLMM) with a single random effect using penalized quasi-likelihood (PQL)
#' 
#' @param formula a two-sided linear formula object describing both the fixed-effects and random-effects part of the model, with the response on the left of a \code{~} operator and the terms, separated by \code{+} operators, on the right.  Random-effects terms are distinguished by vertical bars (\code{|}) separating expressions for design matrices from grouping factors.
#' @param data an optional data frame containing the variables named in
#' @param family a description of the error distribution and link function to be used in the model.  
#' @param weights an optional vector of prior weights with a value for each sample.  
#' @param delta if \code{NULL} estimate delta, if value is given uses this fixed value
#' @param delta.range min and max values (in log space), of the search space for delta to fit the random effect
#' @param maxit max number of PQL iterations
#' @param tol convergence criterion for the 1D search of the delta space
#' @param tol.eta convergence criterion \code{eta} in the PQL iteration
#' @param init.fit \code{fastglmm} object to initialize parameters
#' @param init \code{c("lm", "glm")} method to initialize \code{eta} values
#' @param nthreads number of threads
#'
#' @examples
#' library(MASS)
#' library(lme4)
#' 
#' # GLMM via Laplace approximation
#' fit = glmer(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#' coef(summary(fit))
#' 
#' # GLMM via PQL
#' fit = glmmPQL(y ~ trt + I(week > 2), random = ~ 1 | ID,
#'    family = binomial, data = bacteria, verbose = FALSE)
#' coef(summary(fit))
#' 
#' # GLMM via PQL
#' fit = fastglmm_R(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#' coef(summary(fit))
#
#' @import stats 
#' @importFrom lme4 nobars
#' @importFrom methods is
#' @export
fastglmm_R = function (formula, data, family = gaussian(), weights = NULL, delta = NULL, delta.range = c(-10, 10), maxit = 100, tol = 1e-5, tol.eta = 1e-7, init.fit = NULL, init = c("lm", "glm"), nthreads = 6){

  mc <- match.call()
  init <- match.arg(init)
  
  ## family
  if(is.character(family))
      family <- get(family, mode = "function", envir = parent.frame())
  if(is.function(family)) family <- family()
  if(is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if( !is.null(init.fit) && ! is(init.fit, "fastglmm") ){
    stop("init.fit must be fastglmm object")
  }
  if( maxit < 1 ){
    stop("maxit must be >= 1")
  }

  if( family$family != "poisson" ){
    # glm() handles categorical responses for binomial family
    init <- "glm"
  }

  # decompose formula
  fres = process_formula( formula, data)
  form_fixed = fres$form_fixed
  # update reponse
  form_mod <- update(fres$form_no_offset, zz ~ .)

  if( is.null(weights)) weights <- rep(1, nrow(data))

  # initialize using fixed effects model
  if ( ! is.null(init.fit) ){
    # using previous fastglmm PQL fit
    mf <- model.frame(form_fixed, data)
    offset <- model.offset(mf)
    if( is.null(offset) ) offset <- 0
    # process response thru glm
    suppressWarnings(fit <- glm(form_fixed, family = family, data = data, control=list(maxit=1)))
    # y.orig <- model.response(mf)
    y.orig <- fit$y
    eta.init <- family$linkfun(fitted(init.fit)) + offset

  }else if( init == "glm" ){
    # very slow
    data$weights <- weights
    fit <- glm(form_fixed, family = family, data = data, weights = weights)
    y.orig <- fit$y
    offset <- fit$offset
    if( is.null(offset) ) offset <- 0
    # eta.init = family$linkfun(fitted(fit))
    eta.init = fit$linear.predictors
  }else{
    # faster linear approximation
    # fit lm in transform data
    # PQL part takes longer to converge
    # family$linkfun()
    form <- update(form_fixed, log(.+1e-4) ~ .)
    fit <- lm(form, data = data) # weights
    mf <- model.frame(nobars(formula), data)
    offset <- model.offset(mf)
    y.orig <- model.response(mf)
    eta.init <- fitted(fit)
  }
  # remove names for speed
  names(eta.init) = NULL
  offset = ifelse(is.null(offset), 0, offset)
  eta.init = eta.init - offset

  # print(head(eta.init))

  if( is.null(weights) ){
    w <- weights
  }else{ 
    w <- rep(1, nrow(data))
  }   

  for (i in seq_len(maxit)) {

    # compute linear predictor
    # if change compared to previous is small, break   
    if( i == 1){
      # for glm()
      eta <- eta.init + offset  
    }else{
      # for fastlmm()
      etaold <- eta
      eta <- fitted(fit) + offset  
      if (sum((eta - etaold)^2) < tol.eta){ 
          break          
      }
    }

    # compute updated response and weights 
    mu <- family$linkinv(eta)
    mu.eta.val <- family$mu.eta(eta)
    zz <- eta + (y.orig - mu)/mu.eta.val - offset
    data$zz <- zz
    wz <- w * mu.eta.val^2/family$variance(mu)
    wz <- wz / mean(wz)

    fit <- fastlmm(form_mod, data, 
      weights = wz, 
      delta = delta,
      delta.range = delta.range)  

    fit$formula <- formula 
  }

  fit$response = y.orig
  fit$family <- family
  fit$method <- "PQL"
  fit$iter.pql <- i
  class(fit) <- c("fastglmm", class(fit))

  attr(fit, "call") <- mc
  fit
} 


#' Fit negative binomial mixed model via PQL
#' 
#' Fit negative binomial mixed model (GLMM) with a single random effect using penalized quasi-likelihood (PQL)
#' 
#' @param formula a two-sided linear formula object describing both the fixed-effects and random-effects part of the model, with the response on the left of a \code{~} operator and the terms, separated by \code{+} operators, on the right.  Random-effects terms are distinguished by vertical bars (\code{|}) separating expressions for design matrices from grouping factors.
#' @param data an optional data frame containing the variables named in
#' @param weights an optional vector of prior weights with a value for each sample.  
#' @param maxit max number of NB iterations
#' @param tol convergence criterion for the 1D search of the delta space
#' @param tol.eta convergence criterion \code{eta} in the PQL iteration
#' @param init.fit \code{fastglmm} object to initialize parameters
#' @param init \code{c("lm", "glm")} method to initialize \code{eta} values
#' @param nthreads number of threads
#'
#' @examples
#' library(MASS)
#' library(lme4)
#' 
#' set.seed(101)
#' dd <- expand.grid(f1 = factor(1:3),
#'                f2 = LETTERS[1:2], g=factor(1:9), rep=1:15,
#'        KEEP.OUT.ATTRS=FALSE)
#' mu <- 5*(-4 + with(dd, as.integer(f1) + 4*as.numeric(f2)))
#' dd$y <- rnbinom(nrow(dd), mu = mu, size = 0.5)
#' 
#' # NB GLMM via Laplace approximation
#' fit1 <- glmer.nb(y ~ f1*f2 + (1|g), data=dd)
#' coef(summary(fit1))
#' 
#' # NB GLMM via PQL
#' fit2 <- fastglmm.nb(y ~ f1*f2 + (1|g), data=dd)
#' coef(summary(fit2))
#
#' @importFrom MASS negative.binomial
#' @export
fastglmm_R.nb = function (formula, data, weights, maxit = 100, tol = .Machine$double.eps^0.5, tol.eta = .Machine$double.eps^0.5, init.fit = NULL, init = c("lm", "glm"), nthreads = 6){

  if( missing(weights) ){
    weights = rep(1, nrow(data))
  }

  # fit poisson model
  fit = fastglmm_R(formula, 
        data = data, 
        weights = weights, 
        family = poisson(),
        init.fit = init.fit,
        init = init, 
        nthreads = nthreads, 
        tol = tol)

  # get original counts response
  y.orig = as.numeric(fit$response)

  for(i in seq(maxit)){

    # estimate overdispersion
    theta = nb_theta(y = y.orig, 
              mu = fitted(fit), 
              n = sum(weights),
              weights = weights, 
              left = -5,
              right = 20,
              tol = .Machine$double.eps^0.25)

    ll_prev = logLik(fit)

    # estimate NB model with dispersion fixed
    fit <- fastglmm_R(formula, 
            data = data, 
            weights = weights,
            family = negative.binomial(theta),
            init.fit = fit, 
            nthreads = nthreads, 
            tol = tol, 
            tol.eta = tol.eta)

    # stopping criteria
    if( rel_diff(ll_prev[1], logLik(fit)[1]) < 1e-6) break
  }

  fit$iter.nb = i
  fit
}




rel_diff = function(a,b){
  abs(a - b) / max(abs(a), abs(b))
}
