

library(RUnit)

test_indicator_decomp = function(){
	library(MASS)
	library(fastglmm)
	library(Matrix)
	set.seed(1)

	n = 30000
	ndonors = 20

	info = data.frame(x = rnorm(n))
	info$Indiv = factor(LETTERS[sample(seq(ndonors), n, replace=TRUE)], LETTERS)
	info$Indiv = droplevels(info$Indiv)
	
	library(Matrix)

	# Compute SVD directly from info$Indiv, using count order
	# dsgn = sparse.model.matrix(~ 0 + Indiv, info)
	Z = t(fac2sparse(info$Indiv))
	ZtZ = crossprod(Z)
	dcmp = eigen(ZtZ, symmetric=TRUE)
	# Z = UDV^T
	# Z V^-T D^-1 = U
	# ZtZ = V D U^T U D V^T
	# dcmp2 = svd(Z)

	dcmp$vectors = tcrossprod(Z, solve(dcmp$vectors)) %*% diag(1/sqrt(dcmp$values))
	dcmp$vectors = Matrix(dcmp$vectors, sparse = TRUE)

	dcmp2 = indicator_decomp( info$Indiv )

	checkEqualsNumeric( dcmp$values, dcmp2$values)

	# since there are duplicate eigen-values
	# eigen-vector order can very
	checkEqualsNumeric( dcmp$vectors[,1], dcmp2$vectors[,1])

	C = cor( as.matrix(dcmp$vectors), as.matrix(dcmp2$vectors))
	unique(round(apply(C, 1, max), 3))


	# Weighted rows
	#--------------
	library(Matrix)
	Z = t(fac2sparse(info$Indiv))
	weights = seq(nrow(Z))
	weights = nrow(Z) * (weights / sum(weights))
	ZtZ = crossprod(weights * Z)
	dcmp = eigen(ZtZ, symmetric=TRUE)
	dcmp$vectors = Matrix(zapsmall(dcmp$vectors), sparse = TRUE)
	A = solve(dcmp$vectors)
	D = Diagonal(length(dcmp$values), 1/sqrt(dcmp$values))
	dcmp$vectors = tcrossprod(weights *Z, A) %*% D

	dcmp2 = indicator_decomp( info$Indiv, weights)

	checkEqualsNumeric( dcmp$values, dcmp2$values)
	checkEqualsNumeric( dcmp$vectors[,1], dcmp2$vectors[,1])


}

test_logLik = function(){

	# 1/(sqrt(2 pi) sigma) e^-((x - mu)^2/(2 sigma^2))      

	dmvnorm_r = function(x, mu, Sigma, w=rep(1, ncol(Sigma))){
		r = x - mu
		r = r * sqrt(w) # apply weights
		n = length(r)
		# last term comes from weighting Sigma
		-n/2 * log(2*pi) - 1/2*determinant(Sigma)$modulus[1] - 1/2*crossprod(r, solve(Sigma, r)) - sum(log(1/sqrt(w)))
	}

	library(mvtnorm)
	library(lme4)
	library(RUnit)
	library(fastglmm)

	set.seed(1)
	n = 2000
	y = rnorm(n)
	fit = lm(y ~ 1)
	logLik(fit)

	i = 12
	y_ = y[seq(i)]
	sum(dnorm(y_, 0, sd = 1, log=TRUE))
	dmvnorm(y_, rep(0, i), diag(rep(1,i)), log=TRUE)
	dmvnorm_r(y_, rep(0, i), diag(rep(1,i)))

	i = n
	y_ = y[seq(i)]
	a = mean(y_)
	b = var(y_)
	sum(dnorm(y_, a, sd = sqrt(b), log=TRUE))
	dmvnorm(y_, rep(a, i), diag(rep(b,i)), log=TRUE)
	dmvnorm_r(y_, rep(a, i), diag(rep(b,i)))

	n = 2000
	eta = rnorm(2)
	y = rnorm(n)
	fit = lm(y ~ 1)
	logLik(fit)
	a = mean(y)
	Sigma = diag(rep(var(y)), i)
	dmvnorm(y, rep(a, i), Sigma, log=TRUE)
	dmvnorm_r(y, rep(a, i), Sigma)

	n = 100
	set.seed(1)
	siga = 1e2
	sige = 4
	w = seq(n)
	w = w / mean(w)
	# w[] = 23
	info = data.frame( x = rep(LETTERS[1:20], n/20))
	dsgn = model.matrix(~0+x, info)
	eta = dsgn %*% rnorm(20, sd=sqrt(siga))
	info$y = eta + sqrt(w)*rnorm(n, sd=sqrt(sige))
	info$y = as.numeric(info$y)
	
	fit = lmer(y ~ (1|x), info, weights=w, REML=FALSE)
	logLik(fit)
	summary(fit)
	sigma(fit)^2 / mean(w)^2




	sige_hat = sigma(fit)^2
	siga_hat = attr(VarCorr(fit)[[1]], "stddev")^2
	a = coef(summary(fit))[1,1]
	Sigma = tcrossprod(dsgn)*siga_hat + diag(sige_hat/w)
	dmvnorm(info$y, rep(a, n), Sigma, log=TRUE)
	dmvnorm_r(info$y, rep(a, n), Sigma)


	dcmp = svd(dsgn)
	U = dcmp$u
	d = dcmp$d^2
	# Sigma = U %*% diag(d,2) %*% t(U) + diag(sige_hat/w)
	# Sigma = diag(1/sqrt(w)) %*% (diag(sqrt(w)) %*% U %*% diag(d,2) %*% t(diag(sqrt(w)) %*% U) + diag(rep(sige_hat,n))) %*% diag(1/sqrt(w))
	dmvnorm(info$y, rep(a, n), Sigma, log=TRUE)
	dmvnorm_r(info$y, rep(a, n), Sigma)




	# X = matrix(1, nrow=nrow(info), ncol=1) 
	# Sigma = diag(1/sqrt(w)) %*% (diag(sqrt(w)) %*% U %*% diag(d,2) %*% t(diag(sqrt(w)) %*% U) + diag(rep(sige_hat,n))) %*% diag(1/sqrt(w))
	# r = info$y - X*a

	# crossprod(r, solve(Sigma, r))
	# determinant(Sigma)$modulus[1]


	# Sigma = (diag(sqrt(w)) %*% U %*% diag(d,2) %*% t(diag(sqrt(w)) %*% U) + diag(rep(sige_hat,n))) 
	# crossprod(r*sqrt(w), solve(Sigma, r*sqrt(w)))
	# determinant(diag(1/sqrt(w)) %*% Sigma %*% diag(1/sqrt(w)) )$modulus[1]
	# determinant(Sigma)$modulus[1] + 2*sum(log(1/sqrt(w)))


	# Sigma = (diag(sqrt(w)) %*% U %*% diag(d,2) %*% t(diag(sqrt(w)) %*% U) + diag(rep(sige_hat,n)))
	# Xw = matrix(1, nrow=nrow(info), ncol=1) * sqrt(w)
	# dmvnorm_r(info$y*sqrt(w), Xw*a, Sigma)


	# -n/2 * log(2*pi) - 1/2*determinant(Sigma)$modulus[1] - 1/2*crossprod(r*sqrt(w), solve(Sigma, r*sqrt(w))) - sum(log(1/sqrt(w)))

	# divide X, y by sqrt(w)
	# do SVD of sqrt(w) * 

	sige_hat = sigma(fit)^2
	siga_hat = attr(VarCorr(fit)[[1]], "stddev")^2
	a = coef(summary(fit))[1,1]
	Sigma = tcrossprod(dsgn)*siga_hat + diag(sige_hat/w)
	# dmvnorm_r(yw, rep(a, n), Sigma)




	dcmp = indicator_decomp( factor(info$x), sqrt(w) )
	U = dcmp$vectors
	s = dcmp$values 
	yw = info$y * sqrt(w)
	X = matrix(1, nrow=nrow(info), ncol=1)
	Xw = X * sqrt(w)
	fit2 = fastlmm_R(yw, Xw, U = U, s = s, weights=w)
	fit3 = fastlmm_R(yw, Xw, U = U, s = s, weights=w, delta = fit2$delta)

	fit = lmer(y ~ (1|x), info, weights=w, REML=FALSE)

	logLik(fit)[1]
	fit2$logLik
	fit3$logLik

	cbind(unlist(fit2), unlist(fit3))

	tol = 1e-3
	res = coef(summary(fit))
	checkEqualsNumeric(logLik(fit)[1], fit2$logLik, tol=tol)
	checkEqualsNumeric(sigma(fit)^2, fit2$sig_e, tol=tol)
	checkEqualsNumeric(res[,1], fit2$beta, tol=tol)
	checkEqualsNumeric(res[,2], fit2$beta_se, tol=tol)
	checkEqualsNumeric(VarCorr(fit)[[1]][1], fit2$sig_g, tol=tol)

	VarCorr(fit)[[1]][1]
	fit2$sig_g

	fit2$sig_e

	sige
	siga

	# need to scale residuals, not Y and X
	# r1 =   (X * sqrt(w)) %*% fit2$beta
	# r2 =   X %*% fit2$beta
	# checkEqualsNumeric(r1, r2)


}


test_fastlmm = function(){

	library(MASS)
	library(fastglmm)
	library(Matrix)
	library(RUnit)
	set.seed(1)

	n = 100000
	ndonors = 300

	# n = 3000
	# info = data.frame(x = rnorm(n))
	# info$Indiv = factor(seq(n))

	info = data.frame(x = rnorm(n))
	info$Indiv = factor(sample(seq(ndonors), n, replace=TRUE))
	info$Indiv = droplevels(info$Indiv)
	beta = 1
	eta = 4 + info$x * beta + model.matrix(~ 0 + Indiv, info) %*% rnorm(nlevels(info$Indiv), 0, sqrt(3)) 
	# info$y = rnegbin(n, mu=exp(eta), theta = 10)
	info$y = eta + rnorm(length(eta))

	dcmp = indicator_decomp( info$Indiv )
	X = model.matrix( ~ x, info)

	isSame = function(fit1, fit2){
		ids = intersect(names(fit1), names(fit2))
		ids = ids[-which(ids == 'iter')]
		a = lapply(ids, function(id){
			# cat(id, "...\n")
			checkEqualsNumeric( fit1[[id]], fit2[[id]], tol =  .Machine$double.eps^0.2)
		})
	}

	# devtools::reload("/Users/gabrielhoffman/workspace/repos/fastglmm")

	# 1 response, matrix dcmp$vectors
	U = as.matrix(dcmp$vectors)
	s = dcmp$values
	Y = as.numeric(info$y)
	fit1 = fastlmm_R( Y, X, U = U, s = s)
	fit2 = fastlmm(Y, X, U = U, s = s)
	isSame(fit1, fit2)

	# 1 response, Sparse dcmp$vectors
	U = dcmp$vectors
	s = dcmp$values
	Y = info$y
	fit1 = fastlmm_R( Y, X, U = U, s = s)
	fit2 = fastlmm(Y, X, U = U, s = s)
	isSame(fit1, fit2)

	# multiple responses, matrix dcmp$vectors
	Ym = cbind(Y, runif(length(Y)), runif(length(Y)), runif(length(Y)), runif(length(Y)))
	U = as.matrix(dcmp$vectors)
	s = dcmp$values
	fit1 = fastlmm_R( Ym[,2], X, U = U, s = s)
	fit2 = fastlmm(Ym[,2], X, U = U, s = s)
	isSame(fit1, fit2)

	# multiple responses, Sparse dcmp$vectors
	Ym = cbind(Y, runif(length(Y)), runif(length(Y)), runif(length(Y)), runif(length(Y)), runif(length(Y)), runif(length(Y)), runif(length(Y)), runif(length(Y)))
	U = dcmp$vectors
	s = dcmp$values
	fit1 = fastlmm_R( Ym[,2], X, U = U, s = s)
	fit2 = fastlmm(Ym[,2], X, U = U, s = s)
	isSame(fit1, fit2)

	# batch, matrix dcmp$vectors
	U = as.matrix(dcmp$vectors)
	fitList1 = lapply(seq(ncol(Ym)), function(i){
		fastlmm( Ym[,i], X, U = U, s = s)})
	fitList2 = fastlmm(Ym, X, U = U, s = s)
	res = lapply(seq(ncol(Ym)), function(i){
		isSame(fitList1[[i]], fitList2[[i]])
	})
	checkTrue(unique(unlist(res)))

	# batch, sparse dcmp$vectors
	U = dcmp$vectors
	fitList1 = lapply(seq(ncol(Ym)), function(i){
		fastlmm( Ym[,i], X, U = U, s = s)})
	fitList2 = fastlmm(Ym, X, U = U, s = s)
	res = lapply(seq(ncol(Ym)), function(i){
		isSame(fitList1[[i]], fitList2[[i]])
	})
	checkTrue(unique(unlist(res)))


	# Compare to lme4
	#################

	library(lme4)
	fit = lmer(y ~ x + (1|Indiv), info, REML=FALSE)
	dcmp = indicator_decomp( info$Indiv )
	U = dcmp$vectors
	s = dcmp$values
	fit2 = fastlmm(info$y, X, U = U, s = s)
	
	tol = 1e-3
	res = coef(summary(fit))
	checkEqualsNumeric(logLik(fit)[1], fit2$logLik, tol=tol)
	checkEqualsNumeric(sigma(fit)^2, fit2$sig_e, tol=tol)
	checkEqualsNumeric(res[,1], fit2$beta, tol=tol)
	checkEqualsNumeric(res[,2], fit2$beta_se, tol=tol)
	checkEqualsNumeric(VarCorr(fit)[[1]][1], fit2$sig_g, tol=tol)



	# fastlmm is 100x faster than lmer()
	# system.time(
	# 	replicate(100, lmer(y ~ x + (1|Indiv), info, REML=FALSE)))
	# system.time(
	# 	replicate(100, fastlmm(info$y, X, U = U, s = s)))

	# weights
	#--------

	# info$rnd = sample(LETTERS[1:3], nrow(info), replace=TRUE)
	# info$rnd = factor(info$rnd)
	weights = seq(1, nrow(info))
	weights = nrow(info) * (weights / sum(weights))
	# weights[] = 2
	fit = lmer(y ~ x + (1|Indiv), info, REML=FALSE, weights=weights)
	dcmp = indicator_decomp( info$Indiv, sqrt(weights) )
	U = dcmp$vectors
	s = dcmp$values 
	yw = info$y * sqrt(weights)
	Xw = X * sqrt(weights)
	fit2 = fastlmm_R(yw, Xw, U = U, s = s, weights=weights)
	fit3 = fastlmm(yw, Xw, U = U, s = s, weights=weights)

	a = intersect(names(fit2), names(fit3))
	a = a[a!="iter"]
	a = lapply(a, function(x){
		# message(x)
		checkEqualsNumeric(fit2[[x]], fit3[[x]], tol=1e-3)
	})


	tol = 1e-3
	res = coef(summary(fit))
	checkEqualsNumeric(logLik(fit)[1], fit2$logLik, tol=tol)
	checkEqualsNumeric(sigma(fit)^2/mean(weights), fit2$sig_e, tol=tol)
	checkEqualsNumeric(res[,1], fit2$beta, tol=tol)
	checkEqualsNumeric(res[,2], fit2$beta_se, tol=tol)
	checkEqualsNumeric(VarCorr(fit)[[1]][1], fit2$sig_g, tol=tol)

	# VarCorr(fit)[[1]][1]
	# fit2$sig_g


	# logLik(fit)[1]
	# fit2$logLik


	# sigma(fit)^2 / mean(weights) 
	# fit2$sig_e
	

	# yw = info$y * sqrt(weights)
	# Xw = X * sqrt(weights)
	# it1 = lm.fit(Xw, yw)
	# it2 = lm.wfit(X, info$y, weights)
	# fit3 = lm(y ~ x, info, weights=weights)
	# # plot(residuals(it1), residuals(it2)*sqrt(weights))

	# coef(it1)
	# coef(it2)
	# coef(fit3)

	# sd(it1$residuals)
	# sd(it2$residuals * sqrt(weights))
	# sigma(fit3)

	# fit3 = lm(y ~ x, info, weights=weights)
	# fit4 = lm(yw ~ 0 + Xw)

	# coef(fit3)
	# coef(fit4)

	# logLik(fit3)
	# logLik(fit4)

	# sigma(fit3)
	# sigma(fit4)

	# beta = solve(crossprod(X*sqrt(weights)), crossprod(X* sqrt(weights), info$y * sqrt(weights)))
	# r = info$y * sqrt(weights)  - (X* sqrt(weights)) %*%beta
	# sqrt(crossprod(r)[1] / length(r))

	# r = (info$y  -  X%*%beta)* sqrt(weights)
	# sqrt(crossprod(r)[1] / length(r))




}	







