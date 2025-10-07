


#' Distributional Variance
#'
#' Compute distributional variance from the model fit
#' 
#' @param fit model fit
#' @param method use either the \code{"lognormal"} or \code{"trigamma"} formulas from Nakagawa, et al. (2017)
#' 
#' @details In generalized linear (mixed) models, the link function contributes to the coefficient of determination (Nakagawa, et al., 2012, 2017; McKelvey and Zavoina, 1975).  
#' 
#' 1 - Residuals gives the R2 values from \code{performance::r2_nakagawa(..., approximation="trigamma")}.  Using \code{performance::r2_mckelvey()} use the "lognormal" approximation
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
getDistrVar <- function(fit, fit_null, method = c("trigamma", "lognormal")) {

  method <- match.arg(method)

  # get family identifier
  famID <- getFamilyString(family(fit))

  # remove theta in nb:theta
  famID2 <- gsub("^(.+):.*", "\\1", famID)

  if( isCountModel(fit) ){

    # mean term based on Eqn 5.8 of Nakagawa, et al. 2017.
    if( ! missing(fit_null) ){
      lambda <- getLambdaFromNull( fit_null )
    }else{      
      lambda <- getLambda( fit )
    }
  }

  distVar <- switch( famID2, 
    "gaussian/identity" = {

      w <- get_mean_weights(fit)

      sigma(fit)^2 / w
    },
    "poisson/log" = {
      switch(method,
        "lognormal" = log(1 + 1 / lambda),
        "trigamma" = trigamma(lambda))
    },
    "binomial/logit" = (pi^2) / 3,
    "binomial/probit" = 1,
    "quasipoisson/log" =  {

      omega <- dispersion(fit)

      switch(method,
        "lognormal" = log(1 + omega / lambda),
        "trigamma" = trigamma(lambda / omega))
    },
    "quasibinomial/logit" = stop("Link not supported"),
    "quasibinomial/probit" = stop("Link not supported"),
    "nb" = {

      if( famID == "nb" ){
        theta <- fit_null$theta
      }else{
        theta <- as.numeric(gsub("^(.+):(.*)$", "\\2", famID))
      }

      switch(method,
        "lognormal" = log(1 + 1 / lambda + 1 / theta),
        "trigamma" = trigamma(1/(1/lambda + 1/theta)))
    })

  if (is.null(distVar)) {
    stop("glm family/link not supported: ", famID)
  }

  distVar
}

get_mean_weights = function(fit){

  if( is(fit, "lmerMod") ){
    w = mean(weights(fit))
  }else if( !is.null(fit$prior.weights) ){
    w <- mean(fit$prior.weights)
  }else if( is(fit, "lm") ){
    w <- ifelse(is.null(fit$weights), 1, mean(fit$weights))
  }else if( is(fit, "modelFits") ){
    w <- rep(1, length(sigma(fit)))
  }else{
    w <- 1
  }

  w
}


#' Estimate baseline rate from count model
#' 
#' Estimate baseline rate from count model
#' 
#' @param fit model fit
#' @param ... other args
#'
#' @keywords internal
#' @rdname getLambda
#' @export
setGeneric("getLambda", function(fit,...) {
  standardGeneric("getLambda")
})

#' @rdname getLambda
#' @export
setMethod("getLambda", signature(fit = "fastlmm"), 
  function(fit,...) {

  if( ! '(Intercept)' %in% names(coef(fit)) ){
    stop("Intercept term is required for variance partitioning analysis of a Poisson model")
  }

  # refit model with only intecept term and random effect
  fit_null <- refitModel(fit, interceptOnly=TRUE,...)

  # mean term, including offset
  mu <- coef(fit_null) + mean(fit_null$offset)

  # mean term based on Eqn 5.8 of Nakagawa, et al. 2017.
  # Uses both intercept and variance component
  # follows approach of insight:::.variance_distributional
  lambda <- exp(mu + 0.5*fit_null$sigSq_g)
  as.numeric( lambda )
})


#' @rdname getLambda
#' @export
setMethod("getLambda", signature("glm"), 
  function(fit,...) {

  # refit model with only intercept and offset
  fit_null = update(fit, . ~ 1, 
                offset = fit$offset,
                data = fit$data,
                family = fit$family)

  # mean term, including offset
  mu <- coef(fit_null) + 
          ifelse(is.null(fit_null$offset), 0, mean(fit_null$offset))

  exp(as.numeric(mu))
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



#' Return array with variance component estimates
#' 
#' Return array with variance component estimates
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
#' @importFrom lme4 VarCorr
#' @export
setMethod("varianceTerms", signature("merMod"), 
  function(object,...) {
  sapply(VarCorr(object), function(x){
    attr(x, "stddev")^2
  })
})


#' @rdname varianceTerms
#' @importFrom lme4 findbars
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
#' @param fit model fit  
#' @param ... other arguments, not used here
#' @param distr.method use either the \code{"lognormal"} or \code{"trigamma"} formulas from Nakagawa, et al. (2017)
#' 
#' @details
#' The coefficient of determination (i.e. R^2) is 1 - [Residuals fraction].  This matches \code{performance::r2_nakagawa()} and \code{performance::r2_mckelvey()}, except these use the \code{"lognormal"} method.   
#' 
#' @references
#' Nakagawa, Johnson, Schielzeth. 2017.  The coefficient of determination R2 and intra-class correlation coefficient from generalized linear mixed-effects models revisited and expanded. J. R. Soc. Interface 14: 20170213. \doi{10.1098/rsif.2017.0213}
#'
#' Nakagawa, Shinichi, and Holger Schielzeth. "A general and simple method for obtaining R2 from generalized linear mixed‐effects models." Methods in ecology and evolution 4, no. 2 (2013): 133-142. \doi{10.1111/j.2041-210x.2012.00261.x}
#
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
setGeneric("varpart", function(fit, ..., distr.method = c("trigamma", "lognormal")) {
  standardGeneric("varpart")
})

#' @rdname varpart
#' @export
setMethod("varpart", signature("fastlmm"), 
  function(fit, ..., distr.method = c("trigamma", "lognormal")){
  
  .varpart(fit = fit, distr.method = distr.method)
})

#' @rdname varpart
#' @export
setMethod("varpart", signature("fastglmm"), 
  function(fit, ..., distr.method = c("trigamma", "lognormal")){
  
  .varpart(fit = fit, distr.method = distr.method)
})

#' @rdname varpart
#' @export
setMethod("varpart", signature("glm"), 
  function(fit, ..., distr.method = c("trigamma", "lognormal")){
  
  .varpart(fit = fit, distr.method = distr.method)
})

#' @rdname varpart
#' @export
setMethod("varpart", signature("lm"), 
  function(fit, ..., distr.method = c("trigamma", "lognormal")){
  
  .varpart(fit = fit, distr.method = distr.method)
})

#' @rdname varpart
#' @export
setMethod("varpart", signature("merMod"), 
  function(fit, ..., distr.method = c("trigamma", "lognormal")){
  
  .varpart(fit = fit, distr.method = distr.method)
})

.varpart = function(fit, distr.method = c("trigamma", "lognormal")){

  distr.method <- match.arg(distr.method)
  signal_var <- var(predict(fit)) 
  resid_var <- getDistrVar( fit, method = distr.method ) 
  total_var <- signal_var + resid_var
  eta_var <- colVars(predict(fit, type="terms")) 

  # remove intercept since variance is zero
  eta_var <- eta_var[names(eta_var) != "(Intercept)"]

  # if there are random effects
  eta_var <- c(eta_var, varianceTerms(fit) )

  frac_signal <- as.numeric(signal_var / total_var)

  c(eta_var / sum(eta_var) * frac_signal, 
    Residuals = resid_var / total_var)
}




