


test_predict_fitted = function(){

	library(fastglmm)
	library(MASS)
	library(lme4)
	library(RUnit)

	# use large sample size since GLMM Laplace and PQL
	# converge with large sample size

	set.seed(101)
	dd <- expand.grid(f1 = factor(1:5),
	            f2 = LETTERS[1:5], g=factor(1:30), rep=1:15,
	    KEEP.OUT.ATTRS=FALSE)
	mu <- 50*(-4 + with(dd, as.integer(f1) + 4*as.numeric(f2)))
	dd$y <- rnbinom(nrow(dd), mu = mu, size = 1e8)

	# fastlmm
	##########

	fit2 <- fastlmm(y ~ f1*f2 + (1|g), data=dd)
	fit1 <- lmer(y ~ f1*f2 + (1|g), data=dd, REML=FALSE)

	# coef estimates
	checkEqualsNumeric( coef(summary(fit1))[,1], 
						coef(summary(fit2))[,1], tol=1e-4 )

	# se estimates
	checkEqualsNumeric( coef(summary(fit1))[,2], 
						coef(summary(fit2))[,2], tol=1e-4 )

	checkEqualsNumeric( fitted(fit1), fitted(fit2), tol=1e-4 )
	checkEqualsNumeric( predict(fit1), predict(fit2), tol=1e-4 )
		
	# se estimates
	checkEqualsNumeric( coef(summary(fit1))[,2], 
						coef(summary(fit2))[,2], tol=1e-4 )

	# fastglmm
	##########

	fit2 <- fastglmm(y ~ f1*f2 + (1|g), data=dd, family=poisson())
	fit1 <- glmer(y ~ f1*f2 + (1|g), data=dd, family=poisson())
	fit3 = glmmPQL( y ~ f1*f2 , random = ~ 1|g, dd, family = poisson(), niter=100)

	# coef estimates
	checkEqualsNumeric( coef(summary(fit1))[,1], 
						coef(summary(fit2))[,1], tol=1e-7 )
	checkEqualsNumeric( coef(summary(fit1))[,1], 
						coef(summary(fit3))[,1], tol=1e-7 )

	# se estimates
	checkEqualsNumeric( coef(summary(fit1))[,2], 
						coef(summary(fit2))[,2], tol=1e-2 )
	checkEqualsNumeric( coef(summary(fit1))[,2], 
						coef(summary(fit3))[,2], tol=1e-2 )

	checkEqualsNumeric( fitted(fit1), fitted(fit2), tol=1e-4 )
	checkEqualsNumeric( fitted(fit1), exp(fitted(fit3)), tol=1e-4 )

	checkEqualsNumeric( predict(fit1, type="resp"), 
						predict(fit2, type="resp"), tol=1e-4 )
	checkEqualsNumeric( predict(fit1, type="resp"), 
						predict(fit3, type="resp"), tol=1e-4 )

	checkEqualsNumeric( predict(fit1, type="link"), 
						predict(fit2, type="link"), tol=1e-4 )
	checkEqualsNumeric( predict(fit1, type="link"), 
						predict(fit3, type="link"), tol=1e-4 )

	checkEqualsNumeric( predict(fit1), predict(fit2), tol=1e-5 )
	checkEqualsNumeric( predict(fit1), predict(fit3), tol=1e-5 )

	checkEqualsNumeric( fitted(fit1), fitted(fit2), tol=1e-4 )
	checkEqualsNumeric( suppressWarnings(hatvalues(fit1)), hatvalues(fit2), tol=1e-1 )
	checkEqualsNumeric( unlist(ranef(fit1)), unlist(ranef(fit2)), tol=1e-4 )

	checkEqualsNumeric( vcov(fit1), vcov(fit2), tol=1e-3 )

	checkEqualsNumeric(	residuals(fit1, "response"), 
											residuals(fit2, "response"),
											tol = 1e-3)

	checkEqualsNumeric(	residuals(fit1, "pearson"), 
											residuals(fit2, "pearson"),
											tol = 1e-3)
	
	checkEqualsNumeric(	residuals(fit1, "dev"), 
											residuals(fit2, "dev"),
											tol = 1e-3)


	# fastglmm + offset
	###################

	set.seed(1)
	dd$size = rpois(nrow(dd), 100000)

	fit2 <- fastglmm(y ~ f1*f2 + (1|g) + offset(log(size)), data=dd, family=poisson())
	fit1 <- glmer(y ~ f1*f2 + (1|g)+ offset(log(size)), data=dd, family=poisson())
	fit3 = glmmPQL( y ~ f1*f2 + offset(log(size)), random = ~ 1|g, dd, family = poisson(), niter=100)

	# coef estimates
	checkEqualsNumeric( coef(summary(fit1))[,1], 
						coef(summary(fit2))[,1], tol=1e-7 )
	checkEqualsNumeric( coef(summary(fit1))[,1], 
						coef(summary(fit3))[,1], tol=1e-7 )

	# se estimates
	checkEqualsNumeric( coef(summary(fit1))[,2], 
						coef(summary(fit2))[,2], tol=1e-2 )
	checkEqualsNumeric( coef(summary(fit1))[,2], 
						coef(summary(fit3))[,2], tol=1e-2 )

	checkEqualsNumeric( fitted(fit1), fitted(fit2), tol=1e-4 )
	checkEqualsNumeric( fitted(fit1), exp(fitted(fit3)), tol=1e-4 )

	checkEqualsNumeric( predict(fit1, type="resp"), 
						predict(fit2, type="resp"), tol=1e-4 )
	checkEqualsNumeric( predict(fit1, type="resp"), 
						predict(fit3, type="resp"), tol=1e-4 )

	checkEqualsNumeric( predict(fit1, type="link"), 
						predict(fit2, type="link"), tol=1e-4 )
	checkEqualsNumeric( predict(fit1, type="link"), 
						predict(fit3, type="link"), tol=1e-4 )

	checkEqualsNumeric( predict(fit1), predict(fit2), tol=1e-5 )
	checkEqualsNumeric( predict(fit1), predict(fit3), tol=1e-5 )

	checkEqualsNumeric( fitted(fit1), fitted(fit2), tol=1e-4 )
	checkEqualsNumeric( suppressWarnings(hatvalues(fit1)), hatvalues(fit2), tol=1e-2 )
	checkEqualsNumeric( unlist(ranef(fit1)), unlist(ranef(fit2)), tol=1e-4 )

	checkEqualsNumeric( vcov(fit1), vcov(fit2), tol=1e-3 )

	checkEqualsNumeric(	residuals(fit1, "response"), 
											residuals(fit2, "response"),
											tol = 1e-3)

	checkEqualsNumeric(	residuals(fit1, "pearson"), 
											residuals(fit2, "pearson"),
											tol = 1e-3)
	
	checkEqualsNumeric(	residuals(fit1, "dev"), 
											residuals(fit2, "dev"),
											tol = 1e-3)



	# fastglmm.nb: Negative binomial
	################################

	set.seed(101)
	dd <- expand.grid(f1 = factor(1:2),
	            f2 = LETTERS[1:4], g=factor(1:2), rep=1:15,
	    KEEP.OUT.ATTRS=FALSE)
	mu <- 50*(-4 + with(dd, as.integer(f1) + 4*as.numeric(f2)))
	dd$y <- rnbinom(nrow(dd), mu = mu, size = 1e8)


	dd$y <- rnbinom(nrow(dd), mu = mu, size = 5)

	fit2 <- fastglmm.nb(y ~ f1*f2 + (1|g), data=dd, doCoxReid=FALSE)
	fit1 <- glmer.nb(y ~ f1*f2 + (1|g), data=dd)

	a = gsub("Negative Binomial\\((\\S+)\\)", "\\1", family(fit2)$family)
	b = gsub("Negative Binomial\\((\\S+)\\)", "\\1", family(fit1)$family)
	checkEqualsNumeric(as.numeric(a), as.numeric(b), tol=1e-2)

	fam = negative.binomial(as.numeric(a))
	fit3 = glmmPQL( y ~ f1*f2 , random = ~ 1|g, dd, family = fam, niter=200)

	dispersion(fit2)
	summary(fit1)$dispersion

	# coef estimates
	checkEqualsNumeric( coef(summary(fit1))[,1], 
						coef(summary(fit2))[,1], tol=1e-3 )
	checkEqualsNumeric( coef(summary(fit1))[,1], 
						coef(summary(fit3))[,1], tol=1e-3 )

	# se estimates
	checkEqualsNumeric( coef(summary(fit1))[,2], 
						coef(summary(fit2))[,2], tol=5e-2 )
	checkEqualsNumeric( coef(summary(fit1))[,2], 
						coef(summary(fit3))[,2], tol=5e-2 )

	checkEqualsNumeric( fitted(fit1), fitted(fit2), tol=1e-3 )
	checkEqualsNumeric( fitted(fit1), exp(fitted(fit3)), tol=1e-3 )

	checkEqualsNumeric( predict(fit1, type="resp"), 
						predict(fit2, type="resp"), tol=1e-3 )
	checkEqualsNumeric( predict(fit1, type="resp"), 
						predict(fit3, type="resp"), tol=1e-3 )

	checkEqualsNumeric( predict(fit1, type="link"), 
						predict(fit2, type="link"), tol=1e-4 )
	checkEqualsNumeric( predict(fit1, type="link"), 
						predict(fit3, type="link"), tol=1e-4 )

	checkEqualsNumeric( predict(fit1), predict(fit2), tol=1e-4 )
	checkEqualsNumeric( predict(fit1), predict(fit3), tol=1e-4 )
		

	checkEqualsNumeric( fitted(fit1), fitted(fit2), tol=1e-2 )
	checkEqualsNumeric( suppressWarnings(hatvalues(fit1)), hatvalues(fit2), tol=1e-1 )
	checkEqualsNumeric( unlist(ranef(fit1)), unlist(ranef(fit2)), tol=1e-3 )

	checkEqualsNumeric( vcov(fit1), vcov(fit2), tol=1e-2 )

	checkEqualsNumeric(	residuals(fit1, "response"), 
											residuals(fit2, "response"),
											tol = 1e-3)

	checkEqualsNumeric(	residuals(fit1, "pearson"), 
											residuals(fit2, "pearson"),
											tol = 1e-3)
	
	checkEqualsNumeric(	residuals(fit1, "dev"), 
											residuals(fit2, "dev"),
											tol = 1e-3)

}




test_fastglmm = function(){

	library(MASS)
	library(fastglmm)
	library(Matrix)
	library(lme4)
	library(RUnit)
	set.seed(1)

	n = 10000
	ndonors = 100

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
	indicObj = preprocess_indicator( info$Indiv )
	X = model.matrix( ~ x, info)

	formula = y ~ x + (1|Indiv)
	data = info
	eta = X %*% c(2, .0007) + model.matrix(~Indiv + 0, data) %*% rnorm(nlevels(info$Indiv))
	data$y = rpois(nrow(data), exp(eta))
	family = poisson()

	# source("/Users/gabrielhoffman/workspace/repos/fastglmm/R/fastglmm.R")

	fit0 = fastlmm( y ~ x + (1|Indiv), data)
	coef(summary(fit0))

	w = data$y + 1e-4
	fit00 = fastlmm( log(y + 1e-4) ~ x + (1|Indiv), data, weights=w)
	coef(summary(fit00))


	fit1 = fastglmm( y ~ x + (1|Indiv), data, family = poisson())
	coef(summary(fit1))
	fit1$iter.pql

	# fit1a = fastglmm_R( y ~ x + (1|Indiv), data, family = poisson(), init="glm")
	# coef(summary(fit1a))
	# fit1a$iter.pql

	# fit1b = fastglmm_R( y ~ x + (1|Indiv), data, family = poisson(), init="glm")
	# coef(summary(fit1b))
	# fit1b$iter.pql

	fit2 = glmmPQL( y ~ x, random = ~ 1|Indiv, data, family = poisson(), niter=100)
	coef(summary(fit2))


	fit3 = glmer( y ~ x + (1|Indiv), data, family = poisson())
	coef(summary(fit3))



	# Compare PQL methods
	#---------------------
	fit1 = fastglmm( y ~ (1|Indiv), data, family = poisson())
	fit2 = glmmPQL( y ~ 1, random = ~ 1|Indiv, data, family = poisson(), niter=100)
	fit3 = glmer( y ~ (1|Indiv), data, family = poisson())

	
	# fit1$sigSq_g / (fit1$sigSq_g + fit1$sigSq_e)
	# fit1$sigSq_e / (fit1$sigSq_g + fit1$sigSq_e)
	# calcVarPart(fit3)

	checkEqualsNumeric( unlist(ranef(fit1)), unlist(ranef(fit2)), tol=1e-3)
	checkEqualsNumeric( unlist(ranef(fit1)), unlist(ranef(fit3)), tol=1e-2)


	checkEqualsNumeric( predict(fit1), predict(fit2), tol=1e-5 )
	checkEqualsNumeric( predict(fit1, type="resp"), 
						predict(fit2, type="resp"), 
						tol=1e-6 )

	checkEqualsNumeric( fitted(fit1), exp(fitted(fit2)), tol=1e-6 )

	# checkEqualsNumeric( fixef(fit1), fixef(fit2) )
	checkEqualsNumeric( fixef(fit1), fixef(fit3), tol=1e-3 )
	checkEqualsNumeric( vcov(fit1), vcov(fit2), tol=1e-1  )
	checkEqualsNumeric( coef(summary(fit1))[,1:3], coef(summary(fit2))[,c(1,2,4)], tol=1e-2 )

	fit1$sigSq_g
	fit1$sigSq_e

	fit1 = fastlmm( y ~ x + (1|Indiv), data)
	fit2 = nlme::lme( y ~ x, random = ~ 1|Indiv, data)
	checkEqualsNumeric( sigma(fit1), sigma(fit2), tol=1e-4 )


	# Compare PQL methods: NB fixed
	#------------------------
	fam = negative.binomial(2.3)
	fit1 = fastglmm( y ~ (1|Indiv), data, family = fam)
	fit2 = glmmPQL( y ~ 1, random = ~ 1|Indiv, data, family = fam, niter=100)
	fit3 = glmer( y ~ (1|Indiv), data, family = fam)

	
	fit1$sigSq_g / (fit1$sigSq_g + fit1$sigSq_e)
	fit1$sigSq_e / (fit1$sigSq_g + fit1$sigSq_e)
	# calcVarPart(fit3)

	checkEqualsNumeric( unlist(ranef(fit1)), unlist(ranef(fit2)), tol=1e-5)
	checkEqualsNumeric( unlist(ranef(fit1)), unlist(ranef(fit3)), tol=1e-2)


	checkEqualsNumeric( fitted(fit1), exp(fitted(fit2)), tol=1e-5 )
	# checkEqualsNumeric( fixef(fit1), fixef(fit2) )
	checkEqualsNumeric( fixef(fit1), fixef(fit3), tol=1e-3 )
	checkEqualsNumeric( vcov(fit1), vcov(fit2), tol=1e-3  )
	checkEqualsNumeric( coef(summary(fit1))[,1:3], 
		coef(summary(fit2))[,c(1,2,4)], tol=1e-3)

	fit1$sigSq_g
	fit1$sigSq_e

	fit1 = fastlmm( y ~ x + (1|Indiv), data)
	fit2 = nlme::lme( y ~ x, random = ~ 1|Indiv, data)
	checkEqualsNumeric( sigma(fit1), sigma(fit2), tol=1e-4) 
}


test_theta_ml = function(){

	library(MASS)
	library(RUnit)
	library(fastglmm)

	# no weights
	quine.nb <- glm.nb(Days ~ .^2, data = quine)
	res1 = theta.ml(quine$Days, fitted(quine.nb), limit=200, eps=1e-6)
	res2 = nb_theta(quine$Days, fitted(quine.nb), tol=1e-6)
	checkEqualsNumeric( res1, res2, tol=1e-6)
	
	# with weights
	yeast <- data.frame(cbind(numbers = 0:5, fr = c(213, 128, 37, 18, 3, 1)))
	fit <- glm.nb(numbers ~ 1, weights = fr, data = yeast)

	mu <- fitted(fit)
	n = sum(yeast$fr)
	res1 = theta.ml(yeast$numbers, mu, n, weights = yeast$fr, limit=2000, eps=1e-8)
	res2 = fastglmm:::nb_theta(yeast$numbers, mu, n, yeast$fr, tol=1e-8)
	checkEqualsNumeric( res1, res2, tol=1e-6)
}








