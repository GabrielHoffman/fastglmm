


#' Distributional Variance
#'
#' Compute distributional variance from the model fit
#' 
#' @param fit model fit
#' @param method use either the \code{"trigamma"}, \code{"lognormal"} or  \code{"delta"} formulas from Nakagawa, et al. (2017)
#' @param lambda.method use either \code{"parametric"} or \code{"mean"} method to estimate the mean rate for count models
#'
#' @details In generalized linear (mixed) models, the link function contributes to the coefficient of determination (Nakagawa, et al., 2012, 2017; McKelvey and Zavoina, 1975).  
#' 
#' 1 - Residuals gives the R2 values from \code{performance::r2_nakagawa(..., approximation="trigamma")}.  Using \code{performance::r2_mckelvey()} use the "lognormal" approximation
#'
#' For count models, the distributional variance is a function of the mean count rate. Following Eqn 5.8 of Nakagawa, et al. 2017, this can be estimated using parameters of a model including only intercept and random effect terms.  But this requires refitting the model dropping the rest of the fixed effects.  Instead, computing the mean of the observed counts is a fast approximation.
#' 
#' @references
#' Nakagawa, Johnson, Schielzeth. 2017.  The coefficient of determination R2 and intra-class correlation coefficient from generalized linear mixed-effects models revisited and expanded. J. R. Soc. Interface 14: 20170213. \doi{10.1098/rsif.2017.0213}
#'
#' Nakagawa, Shinichi, and Holger Schielzeth. "A general and simple method for obtaining R2 from generalized linear mixed‐effects models." Methods in ecology and evolution 4, no. 2 (2013): 133-142. \doi{10.1111/j.2041-210x.2012.00261.x}
#'
#' McKelvey, R., Zavoina, W. (1975), "A Statistical Model for the Analysis of Ordinal Level Dependent Variables", Journal of Mathematical Sociology 4, S.103–120.
#
#' @importFrom stats family
#' @keywords internal
#' @export
getDistrVar <- function(fit, fit_null, method = c("trigamma", "lognormal", "delta"), lambda.method = c("parametric", "mean")) {

  method <- match.arg(method)
  lambda.method <- match.arg(lambda.method)

  # get family identifier
  famID <- getFamilyString(family(fit))

  # remove theta in nb:theta
  famID2 <- gsub("^(.+):.*", "\\1", famID)

  if( isCountModel(fit) ){

    # mean term based on Eqn 5.8 of Nakagawa, et al. 2017.
    if( ! missing(fit_null) ){
      lambda <- getLambdaFromNull( fit_null )
    }else{      
      lambda <- getLambda( fit, method = lambda.method )
    }
    lambda = mean(lambda)
  }

  distVar <- switch( famID2, 
    "gaussian/identity" = {

      w <- get_mean_weights(fit)

      sigma(fit)^2 / w
    },
    "poisson/log" = {
      method <- match.arg(method)
      switch(method,
        "delta"     = 1 / lambda,
        "lognormal" = log(1 + 1 / lambda),
        "trigamma"  = trigamma(lambda))
    },
    "binomial/logit"  = (pi^2) / 3,
    "binomial/probit" = 1,
    "quasipoisson/log" = {

      method <- match.arg(method)
      omega <- dispersion(fit)

      switch(method,
        "delta"     = omega / lambda,
        "lognormal" = log(1 + omega / lambda),
        "trigamma"  = trigamma(lambda / omega))
    },
    "quasibinomial/logit" = stop("Link not supported"),
    "quasibinomial/probit" = stop("Link not supported"),
    "nb" = {

      method <- match.arg(method)
      theta <- getTheta( fit )

      switch(method,
        "delta"     = 1/lambda + 1/theta,
        "lognormal" = log(1 + 1 / lambda + 1 / theta),
        "trigamma"  = trigamma(1/(1/lambda + 1/theta)))
    })

  if (is.null(distVar)) {
    stop("glm family/link not supported: ", famID)
  }

  distVar
}

#' Evaluate NB variance when lambda is Inf
#' 
#' Evaluate NB variance when lambda is Inf
#' 
#' @param theta NB overdispersion parameter
#' @param method approximation method
#' 
#' @export
#' @keywords internal
noiseVarNB = function(theta, method = c("trigamma", "lognormal", "delta")){

  method <- match.arg(method)

  switch(method,
    "delta"     = 1/theta,
    "lognormal" = log(1 + 1 / theta),
    "trigamma"  = trigamma(theta))
}


get_mean_weights = function(fit){

  w.mean <- 1

  if( is(fit, "fastglmm") ){
    w.mean <- mean(fit$prior.weights)
  }else if( is(fit, "fastlmm") ){ 
    w.mean <- mean(weights(fit))
  }else if( is(fit, "modelFits") ){
    w.mean <- rep(1, length(sigma(fit)))
  }else {
    w <- weights(fit)
    w.mean <- ifelse(is.null(w), 1, mean(w))
  } 

  w.mean
}


#' Estimate baseline rate from count model
#' 
#' Estimate baseline rate from count model
#' 
#' @param fit model fit
#' @param method use either \code{"parametric"} or \code{"mean"} method to estimate the mean rate for count models
#' @param ... other args
#'
#' For count models, the distributional variance is a function of the mean count rate. Following Eqn 5.8 of Nakagawa, et al. 2017, this can be estimated using parameters of a model including only intercept and random effect terms.  But this requires refitting the model dropping the rest of the fixed effects.  Instead, computing the mean of the observed counts is a fast approximation.
#'
#' @keywords internal
#' @rdname getLambda
#' @export
setGeneric("getLambda", function(fit, method = c("parametric", "mean"),...) {
  standardGeneric("getLambda")
})

#' @rdname getLambda
#' @export
setMethod("getLambda", signature(fit = "fastlmm"), 
  function(fit, method = c("parametric", "mean"),...) {

  method <- match.arg(method)

  if( ! '(Intercept)' %in% names(coef(fit)) ){
    stop("Intercept term is required for variance partitioning analysis of a Poisson model")
  }

  if( method == "parametric"){

    # refit model with only intecept term and random effect
    fit_null <- refitModel(fit, interceptOnly=TRUE,...)

    # mean term, including offset
    eta <- coef(fit_null) + fit_null$offset

    # mean term based on Eqn 5.8 of Nakagawa, et al. 2017.
    # Uses both intercept and variance component
    # follows approach of insight:::.variance_distributional
    lambda <- exp(eta + 0.5*fit_null$sigSq_g)
    lambda <- as.numeric( lambda )
  }else{
    # stop("Not enabled")
    lambda <- mean(fit$response)
  }

  lambda
})


#' @rdname getLambda
#' @export
setMethod("getLambda", signature("glm"), 
  function(fit, method = c("parametric", "mean"),...) {

  method <- match.arg(method)

  if( method == "parametric"){
    # refit model with only intercept and offset
    fit_null <- update(fit, . ~ 1, 
                  offset = fit$offset,
                  data = fit$data,
                  family = fit$family)

    # mean term, including offset
    if( is.null(fit_null$offset) ){
      eta <- coef(fit_null)
    }else{
      eta <- coef(fit_null) + fit_null$offset
    }

    lambda <- exp(as.numeric(eta))
  }else{
    # stop("Not enabled")
    lambda <- mean(fit$y)
  }
  lambda
})

#' @importFrom MASS glm.nb
#' @rdname getLambda
#' @export
setMethod("getLambda", signature("negbin"), 
  function(fit, method = c("parametric", "mean"),...) {

  method <- match.arg(method)

  # glm.nb doesn't support offset as arg, or in update()
  # must be in formula
  # so create new formula and run new glm.nb()
  i <- attr(fit$terms, "response")
  resp <- deparse(attr(fit$terms,"variables")[[i+1]])

  if( method == "parametric"){
    form <- as.formula(paste(resp, "~ 0"))

    if( ! is.null(fit$offset) ){
      form <- update(form, . ~ offset(os)) 
      fit$model$os <- fit$offset
    }else{
      form <- update(form, . ~ 1)
    }

    # refit model with only intercept and offset
    fit_null <- glm.nb(form, 
                  data = fit$model, 
                  etastart = fit$linear.predictors,
                  weights = fit$prior.weights)

    # mean term, including offset
    if( is.null(fit_null$offset) ){
      eta <- coef(fit_null)
    }else{
      eta <- coef(fit_null) + fit_null$offset
    }

    lambda <- exp(as.numeric(eta))
  }else{
    # stop("Not enabled")
    lambda <- mean(fit$y)
  }

  lambda
})





#' Estimate baseline rate from count model
#' 
#' Estimate baseline rate from count model
#' 
#' @param fit_null model fit of null
#' @param ... other args
#'
#' @rdname getLambda
#' @keywords internal
#' @export
setGeneric("getLambdaFromNull", function(fit_null,...) {
  standardGeneric("getLambdaFromNull")
})



#' Array of variance component estimates
#' 
#' Array of variance component estimates
#' 
#' @param object model fit
#' @param ... other args
#'
#' @rdname varianceTerms
#' @export
setGeneric("varianceTerms", function(object,...) {
  standardGeneric("varianceTerms")
})



#' @rdname varianceTerms
#' @export
setMethod("varianceTerms", signature("lm"), 
  function(object,...) {
  NULL
})

#' @rdname varianceTerms
#' @export
setMethod("varianceTerms", signature("glm"), 
  function(object,...) {
  NULL
})

#' @rdname varianceTerms
#' @export
setMethod("varianceTerms", signature("negbin"), 
  function(object,...) {
  NULL
})


#' @rdname varianceTerms
#' @importFrom lme4 VarCorr
#' @export
setMethod("varianceTerms", signature("merMod"), 
  function(object,...) {

  vc <- VarCorr(object)

  res <- lapply(names(vc), function(x){
      v <- attr(vc[[x]], "stddev")^2
      names(v) = paste(x, names(v), sep='.')
      v
    })
  names(res) <- NULL
  res <- unlist(res)
  
  names(res) <- gsub("\\.\\(Intercept\\)", "", names(res))

  res
})


#' @rdname varianceTerms
#' @export
setMethod("varianceTerms", signature("glmmTMB"), 
  function(object,...) {

  vc <- VarCorr(object)$cond

  res <- lapply(names(vc), function(x){
      v <- attr(vc[[x]], "stddev")^2
      names(v) = paste(x, names(v), sep='.')
      v
    })
  names(res) <- NULL
  res <- unlist(res)
  
  names(res) <- gsub("\\.\\(Intercept\\)", "", names(res))

  res
})


#' @rdname varianceTerms
#' @importFrom reformulas findbars
#' @export
setMethod("varianceTerms", signature("fastlmm"), 
  function(object,...) {

  # get name of random effect
  id.ranef <- findbars(formula(object))

  if( length(id.ranef) > 1){
    stop("Only 1 random effect is supported")
  }

  value <- object$sigSq_g
  names(value) <- all.vars(id.ranef[[1]])

  value
})


#' Variance Partitioning Analysis
#' 
#' Compute fraction of variance attributable to each variable in regression model.  Also interpretable as the intra-class correlation after correcting for all other variables in the model.
#' 
#' @param fit regression model fit  
#' @param method select method for count models:  \code{"exact"} or \code{"approximate"} from the current work, or \code{"trigamma"}, \code{"lognormal"} or \code{"delta"} formulas from Nakagawa, et al. (2017)
#' @param pseudocount pseudocount used for \code{"exact"} and \code{"approximate"} methods for count models
#' @param p.tail probability threashold for evaluating expectations for \code{"exact"} methods for count models
#' @param ... other arguments, passed to \code{vpOther()} or \code{vpCounts()}
#'
#' @details
#' For linear model, variance fractions are computed based on the sum of squares explained by each component.  For the linear mixed model, the variance fractions are computed by variance component estimates for random effects and sum of squares for fixed effects.
#'
#' For a generalized linear model, the variance fraction also includes the contribution of the link function so that fractions are reported on the linear (i.e. link) scale rather than the observed (i.e. response) scale. For linear regression with an identity link, fractions are the same on both scales.  But for logit or probit links, the fractions are not well defined on the observed scale due to the transformation imposed by the link function.
#'
#' The variance implied by the link function is the variance of the corresponding distribution (Nakagawa, et al. 2013, 2017)
#'
#' logit -> logistic distribution -> variance is \eqn{\pi^\frac{2}{3}}
#'
#' probit -> standard normal distribution -> variance is 1
#'
#' For count models, Nakagawa, et al. (2013, 2017) propose a large-count approximation.  Instead, we use an exact method described in Hoffman, et al (2026).
#'
#' @references
#' Nakagawa, Johnson, Schielzeth. 2017.  The coefficient of determination R2 and intra-class correlation coefficient from generalized linear mixed-effects models revisited and expanded. J. R. Soc. Interface 14: 20170213. \doi{10.1098/rsif.2017.0213}
#'
#' Nakagawa, and Schielzeth. "A general and simple method for obtaining R2 from generalized linear mixed‐effects models." Methods in ecology and evolution 4, no. 2 (2013): 133-142. \doi{10.1111/j.2041-210x.2012.00261.x}
#'
#' Hoffman, et al. Partitioning gene expression variance using count models. In prep.
#'
#' @examples
#' library(MASS)
#' library(lme4)
#' 
#' fit = fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'   family = binomial(), data = bacteria)
#' 
#' varpart(fit)
#
#' @importFrom matrixStats colVars
#' @rdname varpart
#' @export
setGeneric("varpart", function(fit, method = c("exact", "approximate", "trigamma", "lognormal", "delta"), pseudocount = 1.0, p.tail = 1e-4, ...) {
  standardGeneric("varpart")
})

#' @rdname varpart
#' @export
setMethod("varpart", signature("fastlmm"), 
  function(fit, method = c("exact", "approximate", "trigamma", "lognormal", "delta"), pseudocount = 1.0, p.tail = 1e-4, ...){
  
  .varpart(
    fit = fit, 
    method = method, 
    pseudocount = pseudocount,
    p.tail = p.tail,
    ...)
})

#' @rdname varpart
#' @export
setMethod("varpart", signature("fastglmm"), 
  function(fit, method = c("exact", "approximate", "trigamma", "lognormal", "delta"), pseudocount = 1.0, p.tail = 1e-4, ...){
  
  .varpart(
    fit = fit, 
    method = method, 
    pseudocount = pseudocount,
    p.tail = p.tail,
    ...)
})

#' @rdname varpart
#' @export
setMethod("varpart", signature("glm"), 
  function(fit, method = c("exact", "approximate", "trigamma", "lognormal", "delta"), pseudocount = 1.0, p.tail = 1e-4, ...){
  
  .varpart(
    fit = fit, 
    method = method, 
    pseudocount = pseudocount,
    p.tail = p.tail,
    ...)
})

#' @rdname varpart
#' @export
setMethod("varpart", signature("negbin"), 
  function(fit, method = c("exact", "approximate", "trigamma", "lognormal", "delta"), pseudocount = 1.0, p.tail = 1e-4, ...){
  
  .varpart(
    fit = fit, 
    method = method, 
    pseudocount = pseudocount,
    p.tail = p.tail,
    ...)
})


#' @rdname varpart
#' @export
setMethod("varpart", signature("lm"), 
  function(fit, method = c("exact", "approximate", "trigamma", "lognormal", "delta"), pseudocount = 1.0, p.tail = 1e-4, ...){
  
  .varpart(
    fit = fit, 
    method = method, 
    pseudocount = pseudocount,
    p.tail = p.tail,
    ...)
})

#' @rdname varpart
#' @export
setMethod("varpart", signature("merMod"), 
  function(fit, method = c("exact", "approximate", "trigamma", "lognormal", "delta"), pseudocount = 1.0, p.tail = 1e-4, ...){
  
  .varpart(
    fit = fit, 
    method = method, 
    pseudocount = pseudocount,
    p.tail = p.tail,
    ...)
})

#' @rdname varpart
#' @export
setMethod("varpart", signature("glmmTMB"), 
  function(fit, method = c("exact", "approximate", "trigamma", "lognormal", "delta"), pseudocount = 1.0, p.tail = 1e-4, ...){
  
  .varpart(
    fit = fit, 
    method = method, 
    pseudocount = pseudocount,
    p.tail = p.tail,
    ...)
})



.varpart = function(fit, method = c("exact", "approximate", "trigamma", "lognormal", "delta"), pseudocount = 1.0, p.tail = 1e-4, ...){

  method <- match.arg(method)

  if( isCountModel(fit) && method %in% c("exact", "approximate")){
    vp <- vpCounts(fit, 
          method = method,
          pseudocount = pseudocount, 
          p.tail = p.tail,
          ... )
  }else{
    if( method %in% c("exact", "approximate")){
      method <- "trigamma"
    }
    vp <- vpOther( fit, method = method, ...)
  }

  vp
}


#' @importFrom matrixStats weightedVar
vpOther <- function(fit, method = c("trigamma", "lognormal", "delta"),...){

  method <- match.arg(method)

  # get model weights
  if( is(fit, "fastglmm" ) ){
    w <- fit$prior.weights
  }else{
    w <- c(weights(fit))
  }

  if( !is.null(w) ){
    w <- w / mean(w)
  }

  # signal variancebased on covariates
  #   but not offset
  y.pred <- predict(fit)
  if( !is.null( getOffset(fit) ) ){
    y.pred <- y.pred - getOffset(fit) 
  }

  signal_var <- weightedVar(y.pred, w) 
  distr_var <- getDistrVar( fit, 
                  method = method) 
  total_var <- signal_var + distr_var
  eta_var <- apply(predictTerms(fit), 2, function(x) weightedVar(x, w))

  if( length(eta_var) == 0 ){
    stop("models with no variables not supported")
  }

  # remove intercept since variance is zero
  eta_var <- eta_var[names(eta_var) != "(Intercept)"]

  # if there are random effects
  eta_var <- c(eta_var, varianceTerms(fit) )

  frac_signal <- as.numeric(signal_var / total_var)

  # if family is Negative Binomial
  if( isNB(fit) ){
    resid_var <- noiseVarNB( getTheta(fit), method )
    count_var <- distr_var - resid_var

    res <- c(eta_var / sum(eta_var) * frac_signal,
      CountNoise = count_var / total_var,
      Residuals = resid_var / total_var)

  }else{
    signal <- eta_var / sum(eta_var) * frac_signal

    # if eta_var and frac_signal is effectively zero, 
    # set signal to zero and avoid division by zero
    if( sum(eta_var) < 1e-12 && frac_signal < 1e-12 ){
      signal[] <- 0
    }

    res <- c(signal, 
      Residuals = distr_var / total_var)
  }

  res 
}

# Variance partitioning for count models
#
#' @importFrom matrixStats colVars
vpCounts <- function(fit, method = c("exact", "approximate"),pseudocount = 1.0, p.tail = 1e-4){

  method <- match.arg(method)
 
  mu <- predict(fit, type = "response")
  theta <- getTheta(fit)

  # Approx total noise
  if( is.na(theta) ){
    theta <- Inf
  }

  if( method == "exact"){ 
    # Exact variances
    # integration over NB
    # get mean and variance of log(y + c) given mu, theta
    res <- log_moments_nb_mu(mu, theta, "exact", c = pseudocount, p_tail = p.tail)
    var.signal <- res$var.signal
    var.noise <- res$var.noise
    alpha <- res$alpha
  }else{
    # Approximate variances
    var.poisson <- mu / (mu+pseudocount)^2

    if( ! is.finite(theta) ){
      var.overdisp <- 0
    }else{
      var.overdisp <- (mu^2/theta) / (mu+pseudocount)^2
    }

    var.signal <- var(log(mu + pseudocount) - (mu + mu^2/theta)/(2*(mu+pseudocount)^2))[1] # second order
    # var.signal <- var(log(mu + pseudocount))[1] # first order
    var.noise <- mean(var.poisson + var.overdisp)

    # fraction of variance that is Poisson shot noise
    alpha <- mean(var.poisson / (var.poisson + var.overdisp))
  }

  # scale noise variance by the QL dispersion scale
  var.noise <- var.noise * dispersion(fit)

  # total variance
  var.total <- var.signal + var.noise

  # total signal
  rho2.signal <- var.signal / var.total

  # Exact total noise
  rho2.noise <- var.noise / var.total

  # Approximate fractions
  # fraction of variance for each variable on eta scale
  # Apply these fractions to divide rho2.signal
  eta_var <- colVars(predictTerms(fit))

  if( !is.null(getOffset(fit)) ){
    # account for variance due to offset
    eta_var['offset'] <- var(getOffset(fit))
  }
  eta_var <- eta_var[names(eta_var) != "(Intercept)"]
  eta_var <- c(eta_var, varianceTerms(fit) )
  gamma <- eta_var / sum(eta_var)

  # Variance fractions
  frac <- c(rho2.signal*gamma, 
    CountNoise = rho2.noise*alpha,
    Residuals = rho2.noise*(1-alpha)
    )

  if( !is.null(getOffset(fit)) ){
    # account for variance due to offset
    i <- match("offset", names(frac))
    frac <- frac[-i] / sum(frac[-i])
  }

  frac
}

#' @importFrom lme4 getME
# extract offset from model fit
#' @importFrom lme4 getME
getOffset <- function(fit){

  if( is(fit, "merMod")){
    os <- getME(fit, "offset")
  }else{
    os <- fit[['offset']]
  }

  os
}









