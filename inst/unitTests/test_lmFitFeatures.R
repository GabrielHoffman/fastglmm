



test_lmFitFeatures = function(){

	library(fastlmm)
	library(RUnit)

	# devtools::reload("/Users/gabrielhoffman/workspace/repos/fastlmm")

	n = 30000
	p = 10000
	nc = 8
	set.seed(1)
	y = seq(n)
	X = matrix(rnorm(n*p), n, p)
	colnames(X) = seq(p)
	X_design = matrix(rnorm(n*nc), n,nc)
	w = seq(n)
	w = w / mean(w)
	# w[] = 1

	fitList = lapply(seq(p), function(j){
		lm(y ~ 0 + X_design + X[,j], weights=w)
		})
	beta = do.call(rbind, lapply(fitList, coef))
	se = do.call(rbind, lapply(fitList, function(fit){
		coef(summary(fit))[,2]}))

	fit = fastlmm:::lmFitFeatures(y, X_design, X, colnames(X), w)

	checkEqualsNumeric(beta, fit$coef)
	checkEqualsNumeric(se, fit$se)

	
	fit = fastlmm:::lmFitFeatures(y, X_design, X, colnames(X), w)
	fit2 = fastlmm:::lmFitFeatures_preproj(y, X_design, X, colnames(X), w)


	idx = nc + 1
	checkEqualsNumeric(fit$coef[,idx], fit2$coef)
	checkEqualsNumeric(fit$se[,idx], fit2$se)


	system.time(fit <- fastlmm:::lmFitFeatures(y, X_design, X, colnames(X), w))
	system.time(fit2 <- fastlmm:::lmFitFeatures_preproj(y, X_design, X, colnames(X), w))









	# # Pre-projection with full H
	# H = diag(1,n) - X_design %*% solve(crossprod(X_design)) %*% t(X_design)
	# y_proj = H %*% y
	# X_proj = H %*% X

	# beta = solve(crossprod(X_proj[,1])) %*% crossprod(X_proj[,1], y_proj)

	# beta
	# fit$coef[1,]

	# r_proj = y_proj - X_proj[,1] %*% beta
	# s2 = crossprod(r_proj) / (n-2)
	# sqrt(solve(crossprod(X_proj[,1])) *s2)

	# fit$se[1,]


	# H = diag(1,n) - X_design %*% solve(crossprod(X_design)) %*% t(X_design)
	# b = y - y %*% X_design %*% solve(crossprod(X_design)) %*% t(X_design)

	# checkEqualsNumeric(b, y_proj)



}