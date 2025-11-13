

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
fastglmm.nb = function (formula, data, weights = NULL, maxit = 100, tol = .Machine$double.eps^0.5, tol.eta = .Machine$double.eps^0.5, nthreads = 6){

	fastglmm(formula, 
				data = data, 
				weights = weights,
				maxit = maxit,
				family = negative.binomial(NA),
				tol = tol,
				tol.eta = tol.eta,
				nthreads = nthreads)
}





