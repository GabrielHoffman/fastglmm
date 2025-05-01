

# in Rcpp pass family link functions as Function 
# as in https://github.com/jaredhuling/fastglm/blob/9a04daa4a99761fee4fc87ecdb100a530f96b161/R/fit_glm.R#L199


#' @importFrom stats family gaussian
#' @export
family.fastglmm = function(object,...){
    object$family
}


# Given a formula with response, random effect, and offset
# create form_fixed and form_no_offset
# formula = a ~ (1|b) + (1|c) + offset(v)
#
# process_formula(formula)
# $form_fixed
# a ~ offset(v)
#
# $form_no_offset
# a ~ (1 | b) + (1 | c)
#
#' @importFrom lme4 findbars
process_formula = function(formula, data){

	# 1) extract fixed effects
	form_fixed <- nobars(formula)
	str_rnd_only = paste0('(', as.character(findbars(formula)), ')', collapse=' + ')

	# 2) remove offset
	# adapted from MASS::glmmPQL()
	Terms <- if(missing(data)) terms(form_fixed) else terms(form_fixed, data = data)
	offt <- attr(Terms, "offset")
	offvars <- as.character(attr(Terms, "variables"))[offt + 1L]

	# if offset is the only term, make intercept explicit
    tf <- drop.terms(Terms, offt, keep.response = TRUE)
    tl <- attr(tf,"term.labels")
    if( length(tl) == 0) tl <- "1"
    formula2 = reformulate(tl, response = form_fixed[[2L]],
                        intercept = attr(tf, "intercept"),
                        env = environment(form_fixed))

    # add back random effects
    form_no_offset <- update(formula2, paste('. ~ . +', str_rnd_only))

    list( form_fixed = form_fixed, form_no_offset = form_no_offset)
}


#' Fit generalized linear mixed model via PQL
#' 
#' Fit generalized linear mixed model (GLMM) with a single random effect using penalized quasi-likelihood (PQL)
#' 
#' @param formula a two-sided linear formula object describing both the fixed-effects and random-effects part of the model, with the response on the left of a \code{~} operator and the terms, separated by \code{+} operators, on the right.  Random-effects terms are distinguished by vertical bars (\code{|}) separating expressions for design matrices from grouping factors.
#' @param data an optional data frame containing the variables named in
#' @param family a description of the error distribution and link function to be used in the model.  
#' @param weights an optional vector of prior weights with a value for each sample.  
#' @param delta if \code{NULL} estimate delta, if value is given uses this fixed value
#' @param delta.range min and max values (in log space), of the search space for delta to fit the random effect
#' @param maxit max number of PQL iterations
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
#' # GLMM via Laplace approximation
#' fit = glmer(y ~ trt + I(week > 2) + (1 | ID),
#'				family = binomial(), data = bacteria)
#' coef(summary(fit))
#' 
#' # GLMM via PQL
#' fit = glmmPQL(y ~ trt + I(week > 2), random = ~ 1 | ID,
#' 				family = binomial, data = bacteria, verbose = FALSE)
#' coef(summary(fit))
#' 
#' # GLMM via PQL
#' fit = fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'				family = binomial(), data = bacteria)
#' coef(summary(fit))
#
#' @import stats 
#' @importFrom lme4 nobars
#' @importFrom methods is
#' @export
fastglmm = function (formula, data, family = gaussian(), weights = NULL, delta = NULL,  delta.range = c(-10, 10), maxit = 100, tol = .Machine$double.eps^0.5, tol.eta = .Machine$double.eps^0.5, init.fit = NULL, init = c("lm", "glm"), nthreads = 6){

    mc <- match.call()
    init <- match.arg(init)
	
	## family
    if(is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if(is.function(family)) family <- family()
    if(is.null(family$family)) {
		print(family)
		stop("'family' not recognized")
    }
    if( !is.null(init.fit) && ! is(init.fit, "fastglmm") ){
    	stop("init.fit must be fastglmm object")
    }
    if( maxit < 1 ){
    	stop("maxit must be >= 1")
    }

	if( family$family != "poisson" ){
		# glm() handles categorical responses for binomial family
		init <- "glm"
	}

	# decompose formula
	fres = process_formula( formula, data)
	form_fixed = fres$form_fixed
    # update reponse
    form_mod <- update(fres$form_no_offset, zz ~ .)

	if( is.null(weights)) weights <- rep(1, nrow(data))

	# initialize using fixed effects model
	if ( ! is.null(init.fit) ){
		# using previous fastglmm PQL fit
		mf <- model.frame(form_fixed, data)
    	offset <- model.offset(mf)
    	if( is.null(offset) ) offset <- 0
    	# process response thru glm
    	suppressWarnings(fit <- glm(form_fixed, family = family, data = data, control=list(maxit=1)))
		# y.orig <- model.response(mf)
    	y.orig <- fit$y
		eta.init <- family$linkfun(fitted(init.fit)) + offset

	}else if( init == "glm" ){
		# very slow
		data$weights <- weights
		fit <- glm(form_fixed, family = family, data = data, weights = weights)
    	y.orig <- fit$y
    	offset <- fit$offset
    	if( is.null(offset) ) offset <- 0
		# eta.init = family$linkfun(fitted(fit))
		eta.init = fit$linear.predictors
	}else{
		# faster linear approximation
		# fit lm in transform data
		# PQL part takes longer to converge
		# family$linkfun()
		form <- update(form_fixed, log(.+1e-4) ~ .)
		fit <- lm(form, data = data) # weights
		mf <- model.frame(nobars(formula), data)
    	offset <- model.offset(mf)
		y.orig <- model.response(mf)
		eta.init <- fitted(fit)
	}
	# remove names for speed
	names(eta.init) = NULL
	offset = ifelse(is.null(offset), 0, offset)
	eta.init = eta.init - offset

	# print(head(eta.init))

	if( is.null(weights) ){
		w <- weights
	}else{ 
		w <- rep(1, nrow(data))
	}   

    # strictness tolerance of fastlmm() increases with each PQL iteration
    tol.vary = 10^seq(-1, log10(tol), length.out=round(maxit/10))
    tol.vary = c(tol.vary, rep(tol, maxit - length(tol.vary)))
    tol.vary[] = tol

    for (i in seq_len(maxit)) {

        # compute linear predictor
        # if change compared to previous is small, break   
        if( i == 1){
        	# for glm()
	        eta <- eta.init + offset  
        }else{
          	# for fastlmm()
          	etaold <- eta
	        eta <- fitted(fit) + offset  
	        if (sum((eta - etaold)^2) < tol.eta){ 
	            break          
	        }
	    }

		# compute updated response and weights 
        mu <- family$linkinv(eta)
        mu.eta.val <- family$mu.eta(eta)
        zz <- eta + (y.orig - mu)/mu.eta.val - offset
        data$zz <- zz
        wz <- w * mu.eta.val^2/family$variance(mu)
        wz <- wz / mean(wz)

        # fit model
        if( i == 1 ){
        	# run first time
	        fit <- fastlmm(form_mod, data, 
	        	weights = wz, 
	        	delta = delta,
	        	delta.range = delta.range, 
	        	tol = tol.vary[i])
	    }else{
	        # workhorse after initial fastlmm() fit
	    	fit <- fastlmm.fit(
		    	Y = data$zz, 
		    	X = fit$design, 
		    	Z = fit$Z, 
		    	weights = wz,
	        	delta = delta,
		        delta.range = delta.range, 
		        tol = tol.vary[i])
	    }	   
    }

    fit$response = y.orig
    fit$family <- family
    fit$method <- "PQL"
    fit$iter.pql <- i
    class(fit) <- c("fastglmm", class(fit))

   	attr(fit, "call") <- mc
    fit
} 


