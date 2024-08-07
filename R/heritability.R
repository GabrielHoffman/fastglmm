
# 1) create profile log-likelihood surface from delta surface, and get standard error from hessian of hsq surface at hsq_hat
# 2) fast permutations

#' @importFrom numDeriv hessian
heritability = function(fit, Y, X, U, s, Z, method=c("information", "permutation"), nperms=100){

	method = match.arg(method)

	# estimate of hsq
	hsq_hat = 1 - 1/(1 + 1/fit$delta)
	# with(fit, sig_g / (sig_g + sig_e))

	if( method == "information" ){
		Yu = crossprod(U,Y)
		Xu = crossprod(U,X)
		f = function(hsq){
			delta = 1 / (1/hsq - 1)
			-1*ll_R(delta, Y, X, Yu, Xu, U, s)
		}

		# Fisher information at delta_hat
		infor = hessian(f, hsq_hat)

		# variance of estimate
		se_hsq = sqrt(1 / infor)

		p.value = pnorm(0, hsq_hat, sd=se_hsq, lower.tail=TRUE)

		res = data.frame(hsq = hsq_hat, se_hsq, p.value)

	}else if( method == "permutation"){

		Y_mat = lapply(seq(nperms), function(i){
			sample(Y, length(Y), replace=TRUE)
		})
		Y_mat = do.call(cbind, Y_mat)

		# Run as batch
		fitList = fastlmm.fit(Y_mat, X, Z=Z)

		h_sq_null = sapply(fitList, function(fit){
			1 - 1/(1 + 1/fit$delta)
		})

		# Approximate null distribution with beta
		mu = mean(h_sq_null)
		se = sd(h_sq_null)
		alpha = mu*(mu * (1-mu)/se^2 - 1)
		beta = (1-mu)*(mu * (1-mu)/se^2 - 1)

		# method of moments for beta distribiton
		p.value = pbeta(hsq_hat, alpha, beta, lower.tail=FALSE)

		res = data.frame(hsq = hsq_hat, p.value)
	}

	res
}