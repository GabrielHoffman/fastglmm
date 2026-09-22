
#' Get link function given family
#'
#' Get link function given family
#' 
#' @param family family function
#' 
#' @return link function for model family
#'
#' @examples
#' data(PsychAD)
#'
#' # regression formula
#' form <- PTPRG ~ (1|SubID) + offset(log(libSize))
#'
#' # NB GLMM on PTPRG expression via PQL
#' fit <- fastglmm.nb(form, PsychAD)
#'
#' getLink(family(fit))
#
#' @keywords internal
#' @importFrom stats qnorm
#' @export
getLink = function(family){

  if( !is(family, "family") ){
    stop("Argument must be of class: family")
  }

  famStr <- ifelse(isNB(family), "nb", getFamilyString(family))

  switch( famStr, 
    "gaussian" = identity,
    "gaussian/identity" = identity,
    "poisson/log" = log,
    "quasipoisson/log" = log, 
    "quasibinomial/logit" = function(x) log(x) - log(1-x), 
    "quasibinomial/probit" = function(x) qnorm(x),
    "binomial/logit" = function(x) log(x) - log(1-x),  
    "binomial/probit" = function(x) qnorm(x),
    "nb" = log) 
}

#' Get inverse link function given family
#'
#' Get inverse link function given family
#' 
#' @param family family function
#' 
#' @return inverse link function for model family
#' 
#' @examples
#' data(PsychAD)
#'
#' # regression formula
#' form <- PTPRG ~ (1|SubID) + offset(log(libSize))
#'
#' # NB GLMM on PTPRG expression via PQL
#' fit <- fastglmm.nb(form, PsychAD)
#'
#' getLinkInv(family(fit))
#
#' @keywords internal
#' @importFrom stats pnorm
#' @export
getLinkInv = function(family){

  if( !is(family, "family") ){
    stop("Argument must be of class: family")
  }

  famStr <- ifelse(isNB(family), "nb", getFamilyString(family))

  switch( famStr, 
    "gaussian" = identity,
    "gaussian/identity" = identity,
    "poisson/log" = exp,
    "quasipoisson/log" = exp, 
    "quasibinomial/logit" = function(x) exp(x) / (1+exp(x)),
    "quasibinomial/probit" = function(x) pnorm(x),
    "binomial/logit" = function(x) exp(x) / (1+exp(x)), 
    "binomial/probit" = function(x) pnorm(x),
    "nb" = exp) 
}



