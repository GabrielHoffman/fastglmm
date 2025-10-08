
#' Get link function given family
#'
#' Get link function given family
#' 
#' @param family family function
#' 
#' @keywords internal
#' @importFrom stats qnorm
#' @export
getLink = function(family){
 switch( getFamilyString(family), 
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
#' @keywords internal
#' @importFrom stats pnorm
#' @export
getLinkInv = function(family){

 switch( getFamilyString(family), 
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



