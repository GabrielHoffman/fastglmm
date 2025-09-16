

# Sept 11, 2024
# expand checks of fastglmm.nb

rel_diff = function(a,b){
	abs(a - b) / max(abs(a), abs(b))
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
fastglmm.nb = function (formula, data, weights, maxit = 100, tol = .Machine$double.eps^0.5, tol.eta = .Machine$double.eps^0.5, init.fit = NULL, init = c("lm", "glm"), nthreads = 6){

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




