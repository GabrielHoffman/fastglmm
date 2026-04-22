

# Sept 11, 2024
# expand checks of fastglmm.nb



#' Fit Negative Binomial Mixed Model via PQL
#' 
#' Fit negative binomial mixed model (GLMM) with a single random effect using penalized quasi-likelihood (PQL)
#' 
#' @param formula a two-sided linear formula object describing both the fixed-effects and random-effects part of the model, with the response on the left of a \code{~} operator and the terms, separated by \code{+} operators, on the right.  Random-effects terms are distinguished by vertical bars (\code{|}) separating expressions for design matrices from grouping factors.
#' @param data an optional data frame containing the variables named in
#' @param weights an optional vector of prior weights with a value for each sample.  
#' @param maxit max number of NB iterations
#' @param tol convergence criterion for the 1D search of the delta space
#' @param tol.eta convergence criterion \code{eta} in the PQL iteration
#' @param doCoxReid use Cox-Reid correction for estimating theta in negative binomial model
#' @param nthreads number of threads
#'
#' @examples
#' library(MASS)
#' data(PsychAD)
#'
#' # regression formula
#' form <- PTPRG ~ (1|SubID) + offset(log(libSize))
#'
#' # NB GLMM on PTPRG expression via PQL
#' fit1 <- fastglmm.nb(form, PsychAD)
#' coef(summary(fit1))
#' 
#' # NB GLMM via Laplace approximation
#' # fit2 <- lme4::glmer.nb(form, PsychAD)
#' # coef(summary(fit2))
#
#' @importFrom MASS negative.binomial
#' @export
fastglmm.nb = function (formula, data, weights = NULL, maxit = 100, tol = 1e-3, tol.eta = 1e-3, doCoxReid = nrow(data) < 1000, nthreads = 6){

	fastglmm(formula, 
				data = data, 
				weights = weights,
				maxit = maxit,
				family = negative.binomial(NA),
				tol = tol,
				tol.eta = tol.eta,
				doCoxReid = doCoxReid,
				nthreads = nthreads)
}





