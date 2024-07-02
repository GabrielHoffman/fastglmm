

library(RUnit)

test_indicator_decomp = function(){
	library(MASS)
	library(fastglmm)
	library(Matrix)
	set.seed(1)

	n = 30000
	ndonors = 20

	info = data.frame(x = rnorm(n))
	info$Indiv = factor(sample(seq(ndonors), n, replace=TRUE))
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


}




test_fastlmm = function(){

	library(MASS)
	library(fastglmm)
	library(Matrix)
	library(RUnit)
	set.seed(1)

	n = 100000
	ndonors = 30

	info = data.frame(x = rnorm(n))
	info$Indiv = factor(sample(seq(ndonors), n, replace=TRUE))
	info$Indiv = droplevels(info$Indiv)
	beta = 1
	eta = 4 + info$x * beta + model.matrix(~ 0 + Indiv, info) %*% rnorm(nlevels(info$Indiv), 0, 3) 
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


}	







