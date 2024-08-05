

#' Fast Linear Mixed Model with 1 random effect
#'
#' Fit a linear mixed-effects model with 1 random effect with REML or maximum likelihood using the very fast algorithm and implementation 
#'
#' @param formula a two-sided linear formula object describing both the fixed-effects and random-effects part of the model, with the response on the left of a ‘~’ operator and the terms, separated by ‘+’ operators, on the right.  Random-effects terms are distinguished by vertical bars (‘|’) separating expressions for design matrices from grouping factors.
#' @param data an optional data frame containing the variables named in
#' @param REML logical scalar - Should the estimates be chosen to optimize the REML criterion vs ML?
#' @param weights an optional vector of prior weights with a value for each sample.  When the response has multiple columns, a vector of weight can be reused for each respose, or a matrix the same dimension as the responses matrix can weight each response separately.
#' @param delta  if \code{NULL} estimate delta, if value is given used this fixed values
#' @param delta.range min and max values (in log space), of the search space for delta to fit the random effect
#' @param tol convergence criterion for the 1D search of the delta space
#
#' @examples
#' library(lme4)
#' 
#' fit <- fastlmm(Reaction ~ Days + (1 | Subject), sleepstudy)
#' 
#' fit
#' 
#' summary(fit)
#' 
#' # results are identical for lmer(...,REML=FALSE)
#' fit2 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy, REML=FALSE)
#' coef(summary(fit2))
#' 
#' @details fill in 
#' 
#' @references{
#'   \insertRef{hoffman2013}{fastlmm}
#'
#'   \insertRef{lippert2011}{fastlmm}
#' }
#'
# other args
# verbose = 0L, subset, weights = NULL, na.action, offset, contrasts = NULL
#' @importFrom lme4 findbars nobars
#' @import Rdpack
#' @export
fastlmm = function (formula, data = NULL, REML = TRUE, delta = NULL, weights = NULL, delta.range = c(-10, 10), tol = .Machine$double.eps^0.5){

    mc <- mcout <- match.call()

    # simplest way to extract data
    formula <- as.formula(formula, env =, parent.frame(1L))

    # check that formula has exactly 1 random effect
    fb <- findbars(formula)
    if( length(fb) == 0 ){
    	stop("formula must contain exactly 1 random effect, but none were specified")
    }
    if( length(fb) > 1 ){
    	stop("formula must contain exactly 1 random effect, but ", length(fb), " were specified")
    }

    # check that only 1 random effect variable is used
    vs <- all.vars(fb[[1]])
    if( length(vs) != 1){
    	stop("Only one variable can be used in the random effect")
    }

    # formula with only fixed effects
    form.fixed <- nobars(formula)

    # extract data
    mf <- model.frame( form.fixed, data)
    Y <- model.response( mf )
    X <- model.matrix( mf, data )

    # decomposition of random effect variable
    if( ! is.factor( data[[vs]]) ){
    	stop("Random effect variable must be a factor")
    }

    Z = preprocess_indicator( data[[vs]] )

    # fit model
    fit <- fastlmm.fit(
    	Y = Y, 
    	X = X, 
    	Z = Z, 
    	delta 	= delta,
    	rank 	= ncol(Z), 
    	weights = weights, 
    	tol 	= tol, 
    	delta.range = delta.range)

    # return model fit
   	attr(fit, "call") <- mc
	
    fit
}


