
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
	beta <- solve( QXX, as.matrix(QXY))

	# Eval sig_g
	ru 				<- Yu - Xu %*% beta
	r 				<- Y - X %*% beta
	inv_s_delta_ru 	<- inv_s_delta * ru

	Qrr <- crossprod(ru, inv_s_delta_ru) + (crossprod(r)[1] - crossprod(ru)[1])/ delta
	sig_g <- Qrr[1] / n


	-n/2 * log(2*pi*sig_g) - 1/2 * (sum( log(s + delta ) ) + (n-rank) * log(delta)) - n/2
} 


# browser()
# # hatvalues
# # H_2 = I - V^{-1} + X (X^T V^{-1} X)^{-1} X^T V^{-1}
# Iu = crossprod(U,diag(1,nrow(X)))
# cp_X_low_I <- crossprod(X, diag(1,nrow(X))) - crossprod(Xu, Iu)
# inv_s_delta_Iu <- inv_s_delta * Iu
# QXI <- crossprod(Xu, inv_s_delta_Iu) + cp_X_low_I / delta
# s_expand = c(s, rep(0, nrow(X) - length(s)))
# # H = diag(1, nrow(X)) - diag(1/(s_expand+delta)) + X %*% solve(QXX, QXI) 
# sig_e = delta * sig_g
# V = (diag(1, nrow(X)) * sig_e + tcrossprod(Z) * sig_g)
# H = diag(1, nrow(X)) - solve(V) + solve(V) %*% X %*% solve(t(X) %*% solve(V, X)) %*% t(X) %*% solve(V)
# diag(H)

# methods(class = "fastlmm")

# BLUP
# V = (diag(1, nrow(X)) * sig_e + U %*% diag(s ) * sig_g)
# Z'V^{-1}(y-Xb) * sig_g
# crossprod(Z, solve(V, y - X %*% beta)) * sig_g^2
# crossprod(Z*sig_g, solve(V/sig_g, y - X %*% beta))

	# fitted values
###############

# Z'V^{-1}(y-Xb) 
# U s ()

# Zu =  crossprod(U, Z)
# cp_Z_low_r_low <- crossprod(Z, r) - crossprod(Zu, ru)
# inv_s_delta 	<- 1/(s+delta)
# QZr <- crossprod(Zu, inv_s_delta * ru) + cp_Z_low_r_low / delta
# crossprod(crossprod(U, Z), inv_s_delta * ru)
# crossprod(crossprod(U, U %*% diag(s)), inv_s_delta * ru)
# crossprod(diag(s), inv_s_delta * as.matrix(ru))

# plot(ranef(fit)$Indiv[,1], QZr)

# plot(ranef(fit)$Indiv[,1], as.matrix(fit1$ru))


# Z.mod = Z
# cs = colSums(Z.mod^2)
# idx = order(cs, decreasing = TRUE)
# # idx = seq(length(cs))
# cs = cs[idx]
# Z.mod = Z.mod[, idx]
# vectors = Z.mod %*% Diagonal(ncol(Z.mod), 1/sqrt(cs))

# # diag(cor(as.matrix(vectors), as.matrix(U)))
# diag(cor(as.matrix(vectors), svd(Z)$u))

# diag(cor(as.matrix(Z), as.matrix(with(svd(Z), u %*% diag(d)))))[1:4]

# # have same crossprod
# diag(cor(as.matrix(crossprod(Z)), crossprod(as.matrix(with(svd(Z), u %*% diag(d))))))


# # Z'V^{-1}(y-Xb) * sigSq_g
# V = (diag(1, nrow(X)) * sigSq_e + tcrossprod(Z) * sigSq_g)
# crossprod(Z, solve(V, Y - X %*% beta)) * sigSq_g


#' Fit linear mixed model using SVD of covariance 
#'
#' Fit linear mixed model using SVD of covariance to scale to large sample sizes.
#'
#' @param Y response vector
#' @param X matrix of covariates
#' @param U principal components of covariance matrix
#' @param s eigen values from of covariance matrix
#' @param weights vector weights with value for each sample
#' @param Xu pre-transformed X value
#' @param Yu pre-transformed Y value
#' @param delta ratio of variance components estimated using
#' @param sig_a_fixed if \code{FALSE}, estimate \code{sigSq_a} from data
#' @param rank number of of principal components used 
#' 
#' @details Fit a linear mixed model with a single variance component.
#'
#' @return summary statistics for model fit, and hypothesis testing using X_test_lst, if available
#' 
#' @importFrom stats optimize pnorm sd pbeta
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
		ru 				<<- Yu - Xu %*% beta
		r 				<<- Y - X %*% beta

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
				r 		= r, 
				ru 		= ru,
				pValues	= pValues)
	class(res) <- "fastlmm"
	return(res)
}


