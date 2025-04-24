# # Nov 7, 2024

# test_lmFitResponses = function(){

# 	library(fastlmm)
# 	library(RUnit)

# 	# devtools::reload("/Users/gabrielhoffman/workspace/repos/fastlmm")

# 	n = 14
# 	m = 3
# 	nc = 2
# 	set.seed(1)
# 	Y = matrix(rnorm(n*m), m, n)
# 	colnames(Y) = paste0("s", seq(n))
# 	rownames(Y) = paste0("r", seq(m))

# 	X = matrix(rnorm(n*nc), n,nc)
# 	colnames(X) = paste0("V", seq(nc))
# 	W = matrix(runif(n*m), ,n)	

# 	fitList = lapply(seq(m), function(j){
# 		lm(Y[j,] ~ 0 + X, weights=W[j,])
# 		})
# 	beta = do.call(rbind, lapply(fitList, coef))
# 	sig = sapply(fitList, sigma)
# 	se = do.call(rbind, lapply(fitList, function(fit){
# 		coef(summary(fit))[,2]}))
# 	V = do.call(cbind, lapply(fitList, function(x) c(vcov(x))))
# 	res = do.call(cbind, lapply(fitList, residuals))
# 	hat = do.call(cbind, lapply(fitList, hatvalues))

# 	# matrix
# 	fit = lmFitResponses(Y, X, W, detail=3)

# 	checkEqualsNumeric(beta, fit$coef)
# 	checkEqualsNumeric(se, fit$se)
# 	checkEqualsNumeric(sig^2, fit$sigSq)
# 	checkEqualsNumeric(V, fit$vcov)
# 	checkEqualsNumeric(res, fit$residuals)
# 	checkEqualsNumeric(hat, fit$hatvalues)
	
# 	# devtools::reload("/Users/gabrielhoffman/workspace/repos/fastlmm")
# 	# fit = lmFitResponses(Y, X, colnames(Y), W, detail=1)

# 	# library(Matrix)

# 	# # dgeMatrix not allowed
# 	# fit = lmFitResponses(Matrix(Y), X, W, detail=3)


# }



# test_lmFitFeatures = function(){

# 	library(fastlmm)
# 	library(RUnit)

# 	# devtools::reload("/Users/gabrielhoffman/workspace/repos/fastlmm")

# 	n = 10
# 	p = 3
# 	nc = 2
# 	set.seed(1)
# 	y = seq(n)
# 	X = matrix(rnorm(n*p), n, p)
# 	colnames(X) = seq(p)
# 	X_design = matrix(rnorm(n*nc), n,nc)
# 	w = seq(n)
# 	w = w / mean(w)
# 	# w[] = 1

# 	fitList = lapply(seq(p), function(j){
# 		lm(y ~ 0 + X_design + X[,j], weights=w)
# 		})
# 	beta = do.call(rbind, lapply(fitList, coef))
# 	sig = sapply(fitList, sigma)
# 	se = do.call(rbind, lapply(fitList, function(fit){
# 		coef(summary(fit))[,2]}))
# 	V = do.call(cbind, lapply(fitList, function(x) c(vcov(x))))
# 	res = do.call(cbind, lapply(fitList, residuals))
# 	hat = do.call(cbind, lapply(fitList, hatvalues))

# 	fit = lmFitFeatures(y, X_design, X, w, detail=3, FALSE)

# 	checkEqualsNumeric(beta, fit$coef)
# 	checkEqualsNumeric(se, fit$se)
# 	checkEqualsNumeric(sig^2, fit$sigSq)
# 	checkEqualsNumeric(V, fit$vcovStacked)
# 	checkEqualsNumeric(res, fit$residuals)
# 	checkEqualsNumeric(hat, fit$hatvalues)

# 	fit2 = lmFitFeatures(y, X_design, X, w, detail=2, TRUE)

# 	idx = nc + 1
# 	checkEqualsNumeric(fit$coef[,idx], fit2$coef)
# 	checkEqualsNumeric(fit$se[,idx], fit2$se)
# 	checkEqualsNumeric(fit$sigSq, fit2$sigSq)
# 	checkEqualsNumeric(fit$rdf, fit2$rdf)
# 	checkEqualsNumeric(fit$vcovStacked[idx^2,], fit2$vcovStacked)
# 	checkEqualsNumeric(fit$residuals, fit2$residuals)




# 	# devtools::reload("/Users/gabrielhoffman/workspace/repos/fastlmm")
# 	# fit = fastlmm:::lmFitFeatures(y, X_design, X, colnames(X), w, FALSE)
# 	# fit2 = fastlmm:::lmFitFeatures(y, X_design, X, colnames(X), w, TRUE)
# 	# idx = nc + 1
# 	# checkEqualsNumeric(fit$coef[,idx], fit2$coef)
# 	# checkEqualsNumeric(fit$se[,idx], fit2$se)



# 	# system.time(fit <- lmFitFeatures(y, X_design, X, colnames(X), w, 0, FALSE))
# 	# system.time(fit2 <- lmFitFeatures(y, X_design, X, colnames(X), w, 0, TRUE))

# 	# # Pre-projection with full H
# 	# H = diag(1,n) - X_design %*% solve(crossprod(X_design)) %*% t(X_design)
# 	# y_proj = H %*% y
# 	# X_proj = H %*% X

# 	# beta = solve(crossprod(X_proj[,1])) %*% crossprod(X_proj[,1], y_proj)

# 	# beta
# 	# fit$coef[1,]

# 	# r_proj = y_proj - X_proj[,1] %*% beta
# 	# s2 = crossprod(r_proj) / (n-2)
# 	# sqrt(solve(crossprod(X_proj[,1])) *s2)

# 	# fit$se[1,]


# 	# H = diag(1,n) - X_design %*% solve(crossprod(X_design)) %*% t(X_design)
# 	# b = y - y %*% X_design %*% solve(crossprod(X_design)) %*% t(X_design)

# 	# checkEqualsNumeric(b, y_proj)



# }