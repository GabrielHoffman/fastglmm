


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
getDistrVar <- function(fit, method = c("trigamma", "lognormal")) {

  method <- match.arg(method)

  # get family identifier
  famID <- getFamilyString(family(fit))

  # remove theta in nb:theta
  famID2 <- gsub("^(.+):.*", "\\1", famID)

  if( famID2 %in% c("poisson/log", "quasipoisson/log", "nb")){

    # mean term based on Eqn 5.8 of Nakagawa, et al. 2017.
    lambda <- get_lambda( fit )
  }

  distVar <- switch( famID2, 
    "gaussian/identity" = {

      if( !is.null(fit$prior.weights) ){
        w <- mean(fit$prior.weights)
      }else if( is(fit, "lm") ){
        w <- mean(fit$weights)
      }else{
        w <- 1
      }
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

      omega <- summary(fit)$dispersion

      switch(method,
        "lognormal" = log(1 + omega / lambda),
        "trigamma" = trigamma(lambda / omega))
    },
    "quasibinomial/logit" = stop("Link not supported"),
    "quasibinomial/probit" = stop("Link not supported"),
    "nb" = {
      theta <- as.numeric(gsub("^(.+):(.*)$", "\\2", famID))

      switch(method,
        "lognormal" = log(1 + 1 / lambda + 1 / theta),
        "trigamma" = trigamma(1/(1/lambda + 1/theta)))
    })

  if (is.null(distVar)) {
    stop("glm family/link not supported: ", famID)
  }
  distVar
}

get_lambda <- function( fit,...){
  UseMethod("get_lambda")
}

#' @exportS3Method
get_lambda.fastglmm <- function( fit,...){
  if( ! '(Intercept)' %in% names(coef(fit)) ){
    stop("Intercept term is required for variance partitioning analysis of a Poisson model")
  }

  # refit model with only intecept term and random effect
  fit_null <- refitModel(fit, interceptOnly=TRUE,...)

  # mean term based on Eqn 5.8 of Nakagawa, et al. 2017.
  # Uses both intercept and variance component
  # follows approach of insight:::.variance_distributional
  lambda <- with(fit_null, exp(coefficients + 0.5*sigSq_g))
  as.numeric( lambda )
}

#' @exportS3Method
get_lambda.glm <- function( fit,...){

  # fit intercept-only model
  lambda <- exp(coef(update(fit, . ~ 1, 
    data = fit$data,
    family = family(fit)))[1])
  as.numeric( lambda )
}


#' Return array with variance component estimates
#' 
#' Return array with variance component estimates
#' 
#' @param object model fit
#'
#' @rdname varianceTerms
#' @export
varianceTerms <- function(object) {
  UseMethod("varianceTerms")
}

#' @rdname varianceTerms
#' @export
varianceTerms.lm <- function(object){
  NULL
}

#' @rdname varianceTerms
#' @export
varianceTerms.glm <- function(object){
  NULL
}

#' @rdname varianceTerms
#' @importFrom lme4 VarCorr
#' @export
varianceTerms.merMod <- function(object){
  sapply(VarCorr(object), function(x){
    attr(x, "stddev")^2
  })
}

#' @rdname varianceTerms
#' @importFrom lme4 findbars
#' @export
varianceTerms.fastlmm <- function(object){

  # get name of random effect
  id.ranef <- findbars(formula(object))

  if( length(id.ranef) > 1){
    stop("Only 1 random effect is supported")
  }

  value <- object$sigSq_g
  names(value) <- all.vars(id.ranef[[1]])

  value
}

#' Variance Partitioning Analysis
#' 
#' Compute fraction of variance attributable to each variable in regression model.  Also interpretable as the intra-class correlation after correcting for all other variables in the model. Works for models fit with \code{glm()} and \code{fastglmm()}
#' 
#' @param fit model fit  
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
#' @importFrom matrixStats colVars
#' @export
varpart = function(fit, distr.method = c("trigamma", "lognormal")){

  distr.method <- match.arg(distr.method)

  # # prior weights
  # if( is(fit, "fastglmm") || is.null(weights(fit)) ){
  #   wsqrt <- 1
  # }else{
  #   w <- weights(fit)
  #   w[]  =1
  #   wsqrt <- c(sqrt(w) / mean(sqrt(w)))
  # }
  wsqrt = 1

  signal_var <- var(predict(fit) * wsqrt) 
  resid_var <- getDistrVar( fit, method = distr.method ) 
  total_var <- signal_var + resid_var
  eta_var <- colVars(predict(fit, type="terms") * wsqrt) 

  # remove intercept since variance is zero
  eta_var <- eta_var[names(eta_var) != "(Intercept)"]

  # if there are random effects
  eta_var <- c(eta_var, varianceTerms(fit) )

  frac_signal <- signal_var / total_var

  c(eta_var / sum(eta_var) * frac_signal, 
    Residuals = resid_var / total_var)
}



