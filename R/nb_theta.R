

#' Estimate theta of the Negative Binomial
#' 
#' Given the estimated mean vector, estimate \code{theta} overdispersion parameter of the negative binomial distribution.
#' 
#' @param y vector of observed values from the negative binomial.
#' @param mu estimated mean vector from generalized linear model
#' @param n number of data points (defaults to the sum of \code{weights})
#' @param weights sample weights.  If missing, set to 1
#' @param left left boundary for estimating theta on a log scale.  So -10 corresponds to a minimum theta value of exp(-10)
#' @param right right boundary for estimating theta on a log scale
#' @param tol tolerance of optimization
#'  
#' @details Estimate overdispersion parameter for negative binomial distribution using univariate optimization with Brent's method.  For very large datasets, this can be much faster than the Newton-Raphson method used by \code{MASS::theta.ml()}.  This is faster since the \code{lgamma()} function used in the likelihood is faster than the \code{digamma()} and \code{trigamma()} functions used in the score and information steps.  Also, the univariate optimization is performed on \code{log(theta)}, while the Newton-Raphson approach is performed on the original scale of \code{theta}.
#' 
#' @examples
#' library(MASS)
#' 
#' quine.nb <- glm.nb(Days ~ .^2, data = quine)

#' # estimate from MASS package
#' theta.ml(quine$Days, fitted(quine.nb), limit=200, eps=1e-6)
#' 
#' # faster for large-scale data
#' nb_theta(quine$Days, fitted(quine.nb))
#
#' @seealso \code{MASS::theta.ml()}
#' @export
nb_theta = function(y, mu, n, weights, left=-10, right=20, tol = .Machine$double.eps^0.25){

	if( missing(weights) ){
		weights <- rep(1, length(y))
	}

	if( missing(n) ){
		n <- sum(weights)
	}

	.nb_theta(y, mu, n, weights, left, right, tol)
}


