
#' @importFrom lme4 nobars
#' @export
fastglmm = function (formula, data, family, delta = NULL, weights = NULL, delta.range = c(-10, 10), tol = .Machine$double.eps^0.5, tol.eta = .Machine$double.eps^0.5, nthreads = 6){

	stopifnot(is(family, 'family'))

	if( is.null(weights)) weights <- rep(1, nrow(data))

	# initialize using fixed effects model
	data$weights = weights
	fit <- glm(nobars(formula), family = family, data = data, weights = weights)
    y.orig <- fit$y
    w <- fit$prior.weights
    form.mod <- update(formula, zz ~ .)

    # off <- model.offset(mf)
    off <- 0

    for (i in seq_len(10)) {

        # compute updated response and weights    
        if( i == 1){
        	# for glm()
	        eta <- fit$linear.predictors + off  
        }else{
          	# for fastlmm()
          	etaold <- eta
	        eta <- fitted(fit) + off  
	        if (sum((eta - etaold)^2) < tol.eta){ 
	            break          
	        }
	    }

        mu <- family$linkinv(eta)
        mu.eta.val <- family$mu.eta(eta)
        zz <- eta + (y.orig - mu)/mu.eta.val - off
        data$zz <- as.vector(zz)
        wz <- c(w * mu.eta.val^2/family$variance(mu))
        wx = wz / mean(wz)

        # fit model
        fit <- fastlmm(form.mod, data, weights = wz, delta.range = delta.range, tol = tol)
    }

    fit$family <- family
    fit$method <- "PQL"
    fit$iter.pql <- i
    class(fit) <- c("fastglmm", class(fit))

    fit
} 


