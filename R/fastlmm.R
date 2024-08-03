# See Design docs:
# https://docs.google.com/document/d/1xO0Gdb5mlgP_54uWmpB_iilabDxAf_OQNXfz7J0EWc0/edit

# Option to pass in pre-computed: 
# Xu, Yu, cp_X_low, cp_X_low_Y_low

# Add argument for prior weights?

# Approximate df and rdf
#  just add fixed effect rank to random effect rank
#  assumes they are orthogonal

# eventually, call this _fastglmm
# and call it from within the PQL iteration

# can I estimate theta within PQL?

#' fastglmm
#'
#' Description 
#'
#' @name fastglmm
#' @useDynLib fastglmm 
#' @importFrom Rcpp evalCpp
NULL

#' Fit linear mixed model using SVD of covariance 
#'
#' Fit linear mixed model using SVD of covariance to scale to large sample sizes.
#'
#' @param Y response vector, or matrix with responses as _columns_
#' @param X matrix of covariates
#' @param U principal components of covariance matrix
#' @param s eigen values from of covariance matrix
#' @param delta ratio of variance components estimated using
#' @param rank number of of principal components used 
#' @param weights weight of each sample 
#' @param tol tolerance for estimating delta by Brent's method
#' 
#' @details Fit a linear mixed model with a single variance component.
#'
#' @return summary statistics for model fit, and hypothesis testing using X_test_lst, if available
#' 
#' @export
fastlmm <- function( Y, X, indObj, delta=NULL, sig_a_fixed = FALSE, rank=ncol(U), weights = NULL, tol = .Machine$double.eps^0.25){

	# add data checks here 
	if( !is.matrix(Y) ){
		Y <- as.matrix(Y)
	}

	# apply weights
	dcmp <- indicator_decomp( indObj, weights)
	U <- dcmp$vectors
	s <- dcmp$values 

	if( !is.null(weights) ){
		Y <- Y * sqrt(weights)
		X <- X * sqrt(weights)
		is_weights_one <- FALSE
	}else{
		weights <- matrix(1, nrow(Y), ncol(Y))
		is_weights_one <- TRUE
	}

	# internals assume responses are _rows_
	Y <- t(Y)

	if( ncol(Y) != nrow(X) ){
		stop("dimension of Y and X do not match")
	}

	rank <- min( rank, ncol(U) )

	if( rank < ncol(U)){
		U <- U[,seq_len(rank), drop=FALSE]
		s <- abs(s[seq_len(rank), drop=FALSE])
	}

	# if delta is NULL, estimate its value
	# but setting to -1 for C++ call
	delta <- ifelse( is.null(delta), -1, delta)
	
	if( nrow(Y) == 1){
		if( is(U, "sparseMatrix") ){
			res <- .fastlmm_vms( Y = Y, 
								X = X, 
								U = U, 
								s = s,
								weights = weights, 
								delta = delta, 
								tol = tol)
		}else{
			res = .fastlmm_vmm( Y = Y, 
								X = X, 
								U = U, 
								s = s,
								weights = weights, 
								delta = delta, 
								tol = tol)
		}
		# format results
		res$coefficients <- as.numeric(res$coefficients)
		res$se <- as.numeric(res$se)

		names(res$coefficients) <- colnames(X)
		names(res$se) <- colnames(X)
		rownames(res$vcov) <- colnames(X)
		colnames(res$vcov) <- colnames(X)
		res$rank = ncol(X)

		if( ! is_weights_one ){
			res$weights = weights
		}


		# adapt this to be rdf from H
		res$df.residual <- nrow(X) - ncol(X)

		class(res) <- "fastlmm"
	}else{

		if( is(U, "sparseMatrix") ){
			res <- .fastlmm_mms( Y_all = Y, 
								X = X, 
								U = U, 
								s = s,
								weights = weights, 
								delta = delta, 
								tol = tol)
		}else{
			res <- .fastlmm_mmm( Y_all = Y, 
								X = X, 
								U = U, 
								s = s,
								weights = weights,
								delta = delta, 
								tol = tol)
		}

		# format results for each entry in list here

		class(res) <- "fastlmmList"
	}

	res
}

ll_R <- function( delta, Y, X, Yu, Xu, U, s ){			

	info <- QXX <- sig_g <- NA

	n <- nrow(X)
	rank <- nrow(Xu)

	cp_X_low <- crossprod(X) - crossprod(Xu)
	cp_X_low_Y_low <- crossprod(X, Y) - crossprod(Xu, Yu)

	# Eval Beta
	inv_s_delta 	<- 1/(s+delta)
	inv_s_delta_Yu 	<- inv_s_delta * Yu
	inv_s_delta_Xu 	<- inv_s_delta * Xu

	QXX <- crossprod(Xu, inv_s_delta_Xu) + cp_X_low / delta
	QXY <- crossprod(Xu, inv_s_delta_Yu) + cp_X_low_Y_low / delta
	beta <- solve( QXX, QXY)

	# Eval sig_g
	ru 				<- Yu - Xu %*% beta
	r 				<- Y - X %*% beta
	inv_s_delta_ru 	<- inv_s_delta * ru

	Qrr <- crossprod(ru, inv_s_delta_ru) + (crossprod(r)[1] - crossprod(ru)[1])/ delta
	sig_g <- Qrr[1] / n

	-n/2 * log(2*pi*sig_g) - 1/2 * (sum( log(s + delta ) ) + (n-rank) * log(delta)) - n/2
} 


#' Fit linear mixed model using SVD of covariance 
#'
#' Fit linear mixed model using SVD of covariance to scale to large sample sizes.
#'
#' @param Y response vector
#' @param X matrix of covariates
#' @param U principal components of covariance matrix
#' @param s eigen values from of covariance matrix
#' @param delta ratio of variance components estimated using
#' @param rank number of of principal components used 
#' 
#' @details Fit a linear mixed model with a single variance component.
#'
#' @return summary statistics for model fit, and hypothesis testing using X_test_lst, if available
#' 
#' @importFrom stats optimize
#' @export
fastlmm_R <- function( Y, X, U, s, weights = rep(1, nrow(X)), Xu = NULL, Yu = NULL, delta=NULL, sig_a_fixed = FALSE, rank=ncol(U)){

	rank <- min( rank, ncol(U) )

	if( rank < ncol(U)){
		U <- U[,seq_len(rank), drop=FALSE]
		s <- abs(s[seq_len(rank), drop=FALSE])
	}
	if( is.integer(Y) ){
		Y <- as.numeric(Y)
	}

	if( is.null(Xu) ){
		Xu <- crossprod(U, X)
	}
	if( is.null(Yu) ){
		Yu <- crossprod(U, Y)
	}

	log_interval <- c(10, -10) 

	n <- nrow(Y)

	if( is.null(n) ){
		n <- length(Y)
	}	

	cp_X_low <- crossprod(X) - crossprod(Xu)
	cp_X_low_Y_low <- crossprod(X, Y) - crossprod(Xu, Yu)

	beta <- sigSq_g <- QXX <- 1

	i <- 0
	ll <- function( delta_log ){			
		i <<- i + 1
		delta <- exp(delta_log)

		# Eval Beta
		inv_s_delta 	<- 1/(s+delta)
		inv_s_delta_Yu 	<- inv_s_delta * Yu
		inv_s_delta_Xu 	<- inv_s_delta * Xu

		QXX <<- crossprod(Xu, inv_s_delta_Xu) + cp_X_low / delta
		QXY <- crossprod(Xu, inv_s_delta_Yu) + cp_X_low_Y_low / delta
		beta <<- solve( QXX, as.matrix(QXY))

		# Eval sig_g
		ru 				<- Yu - Xu %*% beta
		r 				<- Y - X %*% beta
		# r <- r * sqrt(weights)

		inv_s_delta_ru 	<- inv_s_delta * ru

		if( sig_a_fixed ){
			sigSq_g <<- 1
		}else{
			Qrr <- crossprod(ru, inv_s_delta_ru) + (crossprod(r)[1] - crossprod(ru)[1])/ delta
			sigSq_g <<- Qrr[1] / n
		}

		-n/2 * log(2*pi*sigSq_g) - 1/2 * (sum( log(s + delta ) ) + (n-rank) * log(delta)) - n/2 + 1/2 * sum(log(weights))
	} 

	if( is.null(delta) ){
		result <- optimize( ll, log_interval, maximum=TRUE)

		# delta <- result$maximum
		delta <- exp(result$maximum)
	}

	# Need to evaluate ll(), so that obj values are evaluated
	log_L <- ll( delta_log = log(delta))

	beta <- array(beta, dimnames=list(rownames(beta)))
	sigSq_e <- delta * sigSq_g

	###################################
	# Hypothesis test using Wald test #
	###################################

	S <- solve( QXX ) * sigSq_g
	beta_se <- sqrt(diag(S))

	pValues <- pnorm( abs(beta), 0, beta_se, lower.tail=FALSE)*2
	
	df <- sum(s[seq_len(rank)]/(s[seq_len(rank)]+delta))
	
	res <- list( logLik 	= log_L, 
				coefficients 	= beta,
				se = beta_se, 
				vcov 	= S,
				delta 	= delta, 
				sigSq_g = sigSq_g, 
				sigSq_e = sigSq_e, 
				iter 	= i,
				df 		= df, 
				pValues	= pValues)
	class(res) <- "fastlmm"
	return(res)
}




# fastlmm = function( Y, X, U, s, delta=NULL, sig_a_fixed = FALSE, rank=ncol(U), n.grid=3){

# 	rank = min( rank, ncol(U) )

# 	if( rank < ncol(U)){
# 		U = U[,seq_len(rank), drop=FALSE]
# 		s = abs(s[seq_len(rank), drop=FALSE])
# 	}

# 	Xu = crossprod(U, X)
# 	Yu = crossprod(U, Y)

# 	log_interval = c(10, -10)

# 	n = nrow(Y)

# 	if( is.null(n) ){
# 		n = length(Y)
# 	}	

# 	cp_X_low = crossprod(X) - crossprod(Xu)
# 	cp_X_low_Y_low = crossprod(X, Y) - crossprod(Xu, Yu)

# 	Q_XX = function(delta, obj){
# 		crossprod(Xu, obj$inv_s_delta_Xu) + cp_X_low / delta
# 	}

# 	Q_Xy = function(delta, obj){
# 		crossprod(Xu, obj$inv_s_delta_Yu) + cp_X_low_Y_low / delta
# 	}

# 	Q_rr = function(delta, obj){
# 		crossprod(obj$ru, obj$inv_s_delta_ru) + (crossprod(obj$r)[1] - crossprod(obj$ru)[1])/ delta
# 	}

# 	obj = list()

# 	ll = function( delta, env=parent.frame() ){			

# 		# Eval Beta
# 		env$obj$inv_s_delta 	= 1/(s+delta)
# 		env$obj$inv_s_delta_Yu 	= env$obj$inv_s_delta * Yu
# 		env$obj$inv_s_delta_Xu 	= env$obj$inv_s_delta * Xu

# 		env$obj$beta = solve( Q_XX(delta, env$obj), Q_Xy(delta, env$obj))

# 		# Eval sig_g
# 		env$obj$ru 				= Yu - Xu %*% env$obj$beta
# 		env$obj$r 				= Y - X %*% env$obj$beta
# 		env$obj$inv_s_delta_ru 	= env$obj$inv_s_delta * env$obj$ru

# 		if( sig_a_fixed ){
# 			env$obj$sig_g  = 1
# 		}else{
# 			env$obj$sig_g = Q_rr(delta, env$obj)[1] / n
# 		}

# 		log_L = -n/2 * log(2*pi*env$obj$sig_g) - 1/2 * (sum( log(s + delta ) ) + (n-rank) * log(delta)) - n/2

# 		return( log_L )
# 	} 

# 	delta_optimal = c()
# 	ll_values = c()

# 	if( is.null(delta) ){
# 		left_value = ll( exp( log_interval[1] ) )
# 		right_value = ll( exp( log_interval[2] ) )
		
# 		delta_grid = exp(seq( log_interval[1], log_interval[2], length.out=n.grid))

# 		for( i in seq_len(length(delta_grid)-1) ){

# 			result = optimize( ll, c( delta_grid[i], delta_grid[i+1]), maximum=TRUE)#, tol=0.00000001)
# 			ll_values[i] = result$objective
# 			delta_optimal[i] = result$maximum

# 			warnings()
# 		} 

# 		result$objective = max( ll_values )
# 		result$maximum = delta_optimal[which.max(ll_values)]

# 		if( left_value > result$objective && left_value > right_value){
# 			delta = exp( log_interval[1] )
# 			log_L = left_value
# 		}else if( right_value > result$objective ){
# 			delta = exp( log_interval[2] )
# 			log_L = right_value
# 		}else{
# 			delta = result$maximum
# 			log_L = result$objective
# 		}
# 	}else{
# 		delta_optimal = append(delta_optimal, delta)
# 		ll_values = append(ll_values, ll(delta))
# 	}

# 	# Need to evaluate ll(), so that obj values are evaluated
# 	log_L = ll(delta)

# 	beta = as.matrix(obj$beta)
# 	sig_g = obj$sig_g
# 	sig_e = delta * sig_g

# 	# get residuals
# 	# alpha is BLUP
# 	# alpha = U %*% diag(s/(s+delta)) %*% crossprod(U, Y - X %*% beta)
# 	# alpha = sweep(U, 2, s/(s+delta), FUN='*') %*% crossprod(U, Y - X %*% beta)
# 	# Y_hat = alpha + X %*% beta
# 	# resid = Y - Y_hat

# 	###################################
# 	# Hypothesis test using Wald test #
# 	###################################

# 	S = solve( Q_XX(delta, obj) ) *sig_g
# 	se_beta = sqrt(diag(S))

# 	pValues = pnorm( abs(beta), 0, se_beta, lower.tail=FALSE)*2
	
# 	df = sum(s[seq_len(rank)]/(s[seq_len(rank)]+delta))
	
# 	return(list( ML 	= log_L, 
# 				delta 	= delta, 
# 				ve 	 	= sig_e, 
# 				vg 		= sig_g, 
# 				df 		= df, 
# 				beta 	= beta, 
# 				se_beta = se_beta,
# 				pValues	= pValues 
# 				# resid 	= resid,
# 				# alpha 	= alpha
# 				))
# }



# #' Fit linear mixed model using SVD of covariance 
# #'
# #' Fit linear mixed model using SVD of covariance to scale to large sample sizes.
# #'
# #' @param Y response vector
# #' @param X matrix of covariates
# #' @param U principal components of covariance matrix
# #' @param s eigen values from of covariance matrix
# #' @param delta ratio of variance components estimated using
# #' @param rank number of of principal components used 
# #' @param W_til components to remove from U.  (Not recommended)
# #' 
# #' @details Fit a linear mixed model with a single variance component.
# #'
# #' @return summary statistics for model fit, and hypothesis testing using X_test_lst, if available
# #' 
# #' @importFrom stats optimize
# #' @export
# fastlmm3 = function( Y, X, U, s, delta=NULL, sig_a_fixed = FALSE, rank=ncol(U), W_til = NULL){

# 	rank = min( rank, ncol(U) )

# 	# U = decomp$vectors[,1:rank]
# 	# s = abs(decomp$values[1:rank])
# 	U = U[,seq_len(rank), drop=FALSE]
# 	s = abs(s[seq_len(rank), drop=FALSE])

# 	if( rank > 1 ){
# 	#	U = rev_matrix( U )
# 	# 	s = rev(abs(decomp$values[1:rank]))
# 	}

# 	Xu = crossprod(U, X)
# 	Yu = crossprod(U, Y)

# 	log_interval = c(10, -10)

# 	n = nrow(Y)

# 	if( is.null(n) ){
# 		n = length(Y)
# 	}	

# 	cp_X_low = crossprod(X) - crossprod(Xu)
# 	cp_X_low_Y_low = crossprod(X, Y) - crossprod(Xu, Yu)

# 	if( ! is.null( W_til) ){
# 		Wu = crossprod(U, W_til)
# 		cp_W_low = crossprod(W_til) - crossprod(Wu)
# 		cp_W_low_X_low = crossprod(W_til, X) - crossprod(Wu, Xu)
# 		cp_W_low_Y_low = crossprod(W_til, Y) - crossprod(Wu, Yu)

# 		cp_W_low[] = 0
# 		cp_W_low_X_low[] = 0
# 		cp_W_low_Y_low[] = 0
		
# 		I_c = diag(1, ncol(W_til))
# 	}

# 	Q_XX = function(delta, obj){
# 		crossprod(Xu, obj$inv_s_delta_Xu) + cp_X_low / delta
# 	}

# 	Q_Xy = function(delta, obj){
# 		crossprod(Xu, obj$inv_s_delta_Yu) + cp_X_low_Y_low / delta
# 	}

# 	Q_rr = function(delta, obj){
# 		crossprod(obj$ru, obj$inv_s_delta_ru) + (crossprod(obj$r)[1] - crossprod(obj$ru)[1])/ delta
# 	}

# 	Q_XW = function(delta, obj){
# 		if( is.null( W_til)) return(0);
# 		t(crossprod(Xu, obj$inv_s_delta_Wu)) + cp_W_low_X_low / delta
# 	}

# 	Q_Wy = function(delta, obj){
# 		if( is.null( W_til)) return(0);
# 		crossprod(Wu, obj$inv_s_delta_Yu) + cp_W_low_Y_low / delta
# 	}

# 	Q_Wr = function(delta, obj){
# 		if( is.null( W_til)) return(0);
# 		cp_W_low_r = crossprod(W_til, obj$r) - crossprod(Wu, obj$ru)

# 		crossprod(Wu, obj$inv_s_delta_ru) + cp_W_low_r / delta
# 	}

# 	Q_WW = function(delta, obj){
# 		if( is.null( W_til)) return(0);
# 		crossprod(Wu, obj$inv_s_delta_Wu) + cp_W_low / delta
# 	}

# 	Omega_XX = function(delta, obj){
# 		if( is.null( W_til) )
# 		Q_XX(delta, obj)
# 		else
# 		Q_XX(delta, obj) + obj$eval_Q_XW_WW %*% obj$eval_Q_XW
# 	}

# 	Omega_Xy = function(delta, obj){
# 		if( is.null( W_til) )
# 		Q_Xy(delta, obj) 
# 		else
# 		Q_Xy(delta, obj) + obj$eval_Q_XW_WW %*% Q_Wy(delta, obj)
# 	}

# 	Omega_rr = function(delta, obj){
# 		if( is.null( W_til) )
# 		Q_rr(delta, obj)
# 		else
# 		Q_rr(delta, obj) + crossprod(obj$eval_Q_Wr, obj$solve_I_Q_WW) %*% obj$eval_Q_Wr 
# 	}

# 	obj = list()

# 	ll = function( delta, env=parent.frame() ){			

# 		# Eval Beta
# 		env$obj$inv_s_delta 	= 1/(s+delta)
# 		env$obj$inv_s_delta_Yu 	= env$obj$inv_s_delta * Yu
# 		env$obj$inv_s_delta_Xu 	= env$obj$inv_s_delta * Xu
		
# 		if( ! is.null( W_til) ){
# 			env$obj$inv_s_delta_Wu 	= env$obj$inv_s_delta * Wu
# 			env$obj$solve_I_Q_WW 	= solve( I_c - Q_WW(delta, env$obj) )
# 			env$obj$eval_Q_XW 		= Q_XW(delta, env$obj)
# 			env$obj$eval_Q_XW_WW = crossprod(env$obj$eval_Q_XW, env$obj$solve_I_Q_WW)
# 		}

# 		env$obj$beta = solve( Omega_XX(delta, env$obj) ) %*% Omega_Xy(delta, env$obj)

# 		# Eval sig_g
# 		env$obj$ru 				= Yu - Xu %*% env$obj$beta
# 		env$obj$r 				= Y - X %*% env$obj$beta
# 		env$obj$inv_s_delta_ru 	= env$obj$inv_s_delta * env$obj$ru
# 		if( ! is.null( W_til) ) env$obj$eval_Q_Wr = Q_Wr(delta, env$obj)

# 		if( sig_a_fixed ){
# 			env$obj$sig_g  = 1
# 		}else{
# 			env$obj$sig_g = Omega_rr(delta, env$obj)[1] / n
# 		}

# 		if( is.null( W_til) )
# 		log_L = -n/2 * log(2*pi*env$obj$sig_g) - 1/2 * (sum( log(s + delta ) ) + (n-rank) * log(delta)) - n/2
# 		else
# 		log_L = -n/2 * log(2*pi*env$obj$sig_g) - 1/2 * (sum( log(s + delta ) ) + (n-rank) * log(delta)) - n/2 - 1/2 * determinant( I_c - Q_WW(delta, env$obj))$modulus[1]

# 		return( log_L )
# 	} 

# 	delta_optimal = c()
# 	ll_values = c()

# 	if( is.null(delta) ){
# 		left_value = ll( exp( log_interval[1] ) )
# 		right_value = ll( exp( log_interval[2] ) )
		
# 		delta_grid = exp(seq( log_interval[1], log_interval[2], length.out=100))

# 		for( i in seq_len(length(delta_grid)-1) ){
# 			 #options(warn = i); message("\n warn =", i, "\n")

# 			result = optimize( ll, c( delta_grid[i], delta_grid[i+1]), maximum=TRUE, tol=0.00000001)
# 			ll_values[i] = result$objective
# 			delta_optimal[i] = result$maximum

# 			warnings()
# 		} 

# 	 	# plot(log(delta_optimal), ll_values)

# 		result$objective = max( ll_values )
# 		result$maximum = delta_optimal[which.max(ll_values)]

# 		if( left_value > result$objective && left_value > right_value){
# 			delta = exp( log_interval[1] )
# 			log_L = left_value
# 		}else if( right_value > result$objective ){
# 			delta = exp( log_interval[2] )
# 			log_L = right_value
# 		}else{
# 			delta = result$maximum
# 			log_L = result$objective
# 		}
# 	}else{
# 		delta_optimal = append(delta_optimal, delta)
# 		ll_values = append(ll_values, ll(delta))
# 	}

# 	# Need to evaluate ll(), so that obj values are evaluated
# 	log_L = ll(delta)

# 	# Evaluate at beta delta
# 	#beta = solve( Omega_XX(delta, obj) ) %*% Omega_Xy(delta, obj)
# 	# sig_g = Omega_rr(delta, obj)[1] / n

# 	beta = obj$beta
# 	sig_g = obj$sig_g
# 	sig_e = delta * sig_g

# 	# get residuals
# 	# alpha is BLUP
# 	# alpha = U %*% diag(s/(s+delta)) %*% crossprod(U, Y - X %*% beta)
# 	alpha = sweep(U, 2, s/(s+delta), FUN='*') %*% crossprod(U, Y - X %*% beta)
# 	Y_hat = alpha + X %*% beta
# 	resid = Y - Y_hat

# 	###################################
# 	# Hypothesis test using Wald test #
# 	###################################

# 	S = solve( Omega_XX(delta, obj) ) *sig_g
# 	se_beta = sqrt(diag(S))

# 	pValues = pnorm( abs(beta), 0, se_beta, lower.tail=FALSE)*2
	
# 	# return log-likelihood surface
# 	###############################

# 	# get degrees of freedom as a function of delta
# 	# df_values = unlist( lapply( delta_optimal, function(delta){ sum(s[seq_len(rank)]/(s[seq_len(rank)]+delta)) } ) )

# 	# surface = as.data.frame( cbind(delta_optimal, df_values, ll_values ) )
# 	# colnames(surface) = c("delta", "df", "log_L")	

# 	# df = sum(s[seq_len(rank)]/(s[seq_len(rank)]+delta))
	
# 	return(list( ML 	= log_L, 
# 				delta 	= delta, 
# 				ve 	 	= sig_e, 
# 				vg 		= sig_g, 
# 				# df 		= df, 
# 				beta 	= beta, 
# 				se_beta = se_beta,
# 				pValues	= pValues, 
# 				# surface = surface,
# 				resid 	= resid,
# 				alpha 	= alpha))
# }


