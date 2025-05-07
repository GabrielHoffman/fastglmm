


# cd /sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonsubtype_level/dynamic_state_eQTL/Aging/run_dynamic_eQTL_with_fixed_theta/EN/another/for_expanded2_test/chr2/test

# /sc/arion/projects/CommonMind/zengb02/single_cell_eQTL_for_Gabriel_02212024/eQTL_detection_on_pearsonsubtype_level/dynamic_state_eQTL/Aging/run_dynamic_eQTL_with_fixed_theta/EN/another/for_expanded2_test/chr2/test


# , containing merge phe+geno in “adjust_merged_pheno_geno_NGEF”, covariate infor in “metadata_NGEF”, and pseudotime infor in “pseudo_NGEF”. You can also find the R code I used to run nb test in nb_univariate_Aging_v4_part8.R


# library(tidyverse)
# df = read_tsv("metadata_NGEF") %>%
# 			as.data.frame

# library(fastglmm)
# library(fastlmm)
# library(lme4)
# library(MASS)

# df$y = rpois(nrow(df), 43)
# df$SubID = factor(df$SubID)

# form = y ~ offset(log(nCount_RNA)) + Sex + (1|SubID)
# fit0 = glmer( form, df, family=poisson())
# coef(summary(fit0))



# form = y ~ offset(log(nCount_RNA)) + Sex
# fit1 = glmmPQL( form, random = ~ 1|SubID, df, family=poisson())
# coef(summary(fit1))




# form = y ~ offset(log(nCount_RNA)) + Sex + (1|SubID)
# fit1 = fastglmm( form, df, family=poisson())
# coef(summary(fit1))


# devtools::reload("/Users/gabrielhoffman/workspace/repos/fastglmm")

# form = y ~ Sex + (1|SubID)
# fit1 = fastglmm( form, df, family=poisson(), init="glm")







# form = y ~ Sex + (1|SubID)
# coef(fastlmm( form, df))

# form = y ~ offset(log(nCount_RNA)) + Sex + (1|SubID)
# coef(fastlmm( form, df))


# fixef(lmer(form, df, REML=FALSE))


# form = y ~ Sex + (1|SubID)
# fit1 = fastlmm( form, df)
# coef(summary(fit1))










# library(fastlmm)
# library(fastglmm)
# library(MASS)
# library(lme4)


# info = bacteria
# info$off = seq(nrow(info))/8



# fit = glmmPQL(y ~ trt + offset(off), random = ~ 1 | ID,
#                              family = binomial, data = info)
# coef(summary(fit))



# fit = fastglmm(y ~ trt + (off) + (1 | ID),
#                      family = binomial, data = info, init = "glm")
# head(fitted(fit))



# fit = glmer(y ~ trt + (off) + (1 | ID),
#                      family = binomial, data = info)
# head(fitted(fit))


# # devtools::reload("/Users/gabrielhoffman/workspace/repos/fastglmm")


# form = y ~ trt  + (1 | ID) + offset(1*off)

# fit = fastglmm(form,
#                      family = binomial, data = info, init = "lm")
# coef(summary(fit))


# fit = fastglmm(form,
#                      family = binomial, data = info, init = "glm")
# coef(summary(fit))


# fit = fastglmm(form,
#                      family = binomial, data = info, init.fit=fit)
# coef(summary(fit))






# fit = glmmPQL(nobars(form), random = ~ 1 | ID,
#                              family = binomial, data = info)
# coef(summary(fit))






# library(lme4)
# library(fastlmm)
# library(fastglmm)

# form = Reaction ~ offset(Days) + (1 | Subject)

# fit0 <- lmer(form, sleepstudy, REML=FALSE)
# coef(summary(fit0))

# fit1 <- fastlmm(form, sleepstudy)
# coef(summary(fit1))

# fit2 <- fastglmm(form, sleepstudy, family=gaussian())
# coef(summary(fit2))




# form = Reaction ~ I(Days/10) + offset(Days) + (1 | Subject)

# fam = gaussian(link='identity')

# fit0 <- glmer(form, sleepstudy, family=fam)
# coef(summary(fit0))


# fit1 <- glmmPQL(nobars(form), random = ~1 | Subject, sleepstudy, family=fam)
# coef(summary(fit1))



# fit2 <- fastglmm(form, sleepstudy, family=fam)
# coef(summary(fit2))











# fastglmm.nb.glmer = function (formula, data, weights = NULL, delta = NULL,  delta.range = c(-10, 10), maxit = 100, tol = .Machine$double.eps^0.5, tol.eta = .Machine$double.eps^0.5, init.fit = NULL, init = c("lm", "glm"), nthreads = 6){

# 	# fit poisson model
# 	fit = glmer(formula, data, weights = weights, family = poisson())

# 	for(i in seq(1)){

# 		# estimate overdispersion
# 		theta = theta.ml(y = fit@frame$y, 
# 			 	mu = fitted(fit),
# 			 	weights = fit@resp$weights, limit=20, eps = .Machine$double.eps^0.25)

# 		# estimate NB model with dispersion fixed
# 		fit = glmer(formula, data, weights = weights, family = negative.binomial(theta))

# 		# glmer code to set theta
# 		# https://github.com/lme4/lme4/blob/bfd7a44d0a718fff090412871504858559a0829f/R/nbinom.R#L17C1-L29C2
# 	}

# 	# models uncertainty in theta estimate
# 	# jointly optimizes theta and other paramers
# 	lme4:::optTheta(fit, log(theta) + c(-3,3), tol=5e-05)
# }


# fit1 = fastglmm.nb( y ~ x + (1|Indiv), data)
# fit2 = glmer.nb( y ~ x + (1|Indiv), data)

# fit.tmp = glmer( y ~ x + (1|Indiv), data, family=poisson())
# fit3 = lme4:::optTheta(fit.tmp, c(-3,8), tol=5e-05)


# family(fit1)
# family(fit2)
# family(fit3)

# vcov(fit1)
# vcov(fit2)
# vcov(fit3)

# coef(summary(fit1))
# coef(summary(fit2))
# coef(summary(fit3))

# theta = 1e-5
# a = lme4:::setNBdisp(fit1, theta)
# a = refit(a)
# b = fit1
# b@call$family$theta = theta
# b = refit(b)

# residuals(fit1)[1:3]
# residuals(a)[1:3]



# coef(summary(a))
# coef(summary(b))


# library(glmmTMB)
# library(NBZIMM)



# system.time(replicate(10, fastglmm.nb( y ~ x + (1|Indiv), data)))
# system.time(replicate(10, a<-glmer( y ~ x + (1|Indiv), data, family=poisson())))
# system.time(replicate(10, glmer.nb( y ~ x + (1|Indiv), data)))
# system.time(replicate(10, glmmTMB( y ~ x + (1|Indiv), data, family=nbinom2)))
# system.time(replicate(10, glmm.nb(y ~ x, random = ~ 1|Indiv, data, verbose=FALSE)))




# fit2 = glmer.nb( y ~ x + (1|Indiv), data)



# fam = negative.binomial(1e5)
# fit1 = fastglmm( y ~ x + (1|Indiv), data, family=fam)
# coef(summary(fit1))

# fam = negative.binomial(1e-5)
# fit1 = fastglmm( y ~ x + (1|Indiv), data, family=fam)
# coef(summary(fit1))


# fam = negative.binomial(1e5)
# fit2 = glmer( y ~ x + (1|Indiv), data, family=fam)
# coef(summary(fit2))

# fam = negative.binomial(1e-5)
# fit2 = glmer( y ~ x + (1|Indiv), data, family=fam)
# coef(summary(fit2))



# fit2 = glmer.nb( y ~ x + (1|Indiv), data)
# coef(summary(fit2))

# fit3 = glmer( y ~ x + (1|Indiv), data, family=family(fit2))
# coef(summary(fit3))



# library(MASS)
# library(fastglmm)
# library(lme4)
# library(glmmTMB)


# set.seed(101)
# dd <- expand.grid(f1 = factor(1:3),
#                f2 = LETTERS[1:2], g=1:20, rep=1:15,
#        KEEP.OUT.ATTRS=FALSE)
# mu <- 5*(-4 + with(dd, as.integer(f1) + 4*as.numeric(f2)))
# dd$y <- rnbinom(nrow(dd), mu = mu, size =1)
# dd$g = factor(dd$g)
# nrow(dd)

# fit1 <- fastglmm.nb(y ~ f1*f2 + (1|g), data=dd)
# fit2 <- glmer.nb(y ~ f1*f2 + (1|g), data=dd)
# fit3 <- fastglmm.nb.glmer(y ~ f1*f2 + (1|g), data=dd)
# fit4 <- glmer(y ~ f1*f2 + (1|g), data=dd, family=family(fit2))
# fit5 <- glmmTMB(y ~ f1*f2 + (1|g), data=dd, family=nbinom2)


# family(fit1)
# family(fit2)
# family(fit3)
# family(fit4)
# family(fit5)

# diag(vcov(fit1))
# diag(vcov(fit2))
# diag(vcov(fit3))
# diag(vcov(fit4))

# coef(summary(fit1))
# coef(summary(fit2))
# coef(summary(fit3))
# coef(summary(fit4))
# coef(summary(fit5))$cond



# fitted(fit1)[1:3]
# fitted(fit2)[1:3]
# fitted(fit3)[1:3]
# fitted(fit4)[1:3]


# family(fastglmm.nb(y ~ f1*f2 + (1|g), data=dd))
# family(fastglmm.nb(y ~ f1*f2 + (1|g), data=dd, maxit=2))


test_predict_fitted = function(){

	library(fastlmm)
	library(MASS)
	library(lme4)
	library(RUnit)

	set.seed(101)
	dd <- expand.grid(f1 = factor(1:3),
	            f2 = LETTERS[1:2], g=factor(1:9), rep=1:15,
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
						coef(summary(fit3))[,2], tol=1e-3 )

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
		

	# fastglmm
	##########

	dd$y <- rnbinom(nrow(dd), mu = mu, size = 5)

	fit2 <- fastglmm.nb(y ~ f1*f2 + (1|g), data=dd)
	fit1 <- glmer.nb(y ~ f1*f2 + (1|g), data=dd)

	a = gsub("Negative Binomial\\((\\S+)\\)", "\\1", family(fit2)$family)
	b = gsub("Negative Binomial\\((\\S+)\\)", "\\1", family(fit1)$family)
	checkEqualsNumeric(as.numeric(a), as.numeric(b), tol=1e-2)

	fam = negative.binomial(as.numeric(a))
	fit3 = glmmPQL( y ~ f1*f2 , random = ~ 1|g, dd, family = fam, niter=200)

	# coef estimates
	checkEqualsNumeric( coef(summary(fit1))[,1], 
						coef(summary(fit2))[,1], tol=1e-3 )
	checkEqualsNumeric( coef(summary(fit1))[,1], 
						coef(summary(fit3))[,1], tol=1e-3 )

	# se estimates
	checkEqualsNumeric( coef(summary(fit1))[,2], 
						coef(summary(fit2))[,2], tol=1e-2 )
	checkEqualsNumeric( coef(summary(fit1))[,2], 
						coef(summary(fit3))[,2], tol=1e-2 )

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
		

}




test_fastglmm = function(){

	library(MASS)
	library(fastlmm)
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

	fit1a = fastglmm_R( y ~ x + (1|Indiv), data, family = poisson(), init="glm")
	coef(summary(fit1a))
	fit1a$iter.pql

	fit1b = fastglmm_R( y ~ x + (1|Indiv), data, family = poisson(), init="glm")
	coef(summary(fit1b))
	fit1b$iter.pql

	fit2 = glmmPQL( y ~ x, random = ~ 1|Indiv, data, family = poisson(), niter=100)
	coef(summary(fit2))


	fit3 = glmer( y ~ x + (1|Indiv), data, family = poisson())
	coef(summary(fit3))



	# Compare PQL methods
	#---------------------
	fit1 = fastglmm( y ~ (1|Indiv), data, family = poisson())
	fit2 = glmmPQL( y ~ 1, random = ~ 1|Indiv, data, family = poisson(), niter=100)
	fit3 = glmer( y ~ (1|Indiv), data, family = poisson())

	
	fit1$sigSq_g / (fit1$sigSq_g + fit1$sigSq_e)
	fit1$sigSq_e / (fit1$sigSq_g + fit1$sigSq_e)
	# calcVarPart(fit3)

	checkEqualsNumeric( ranef(fit1), unlist(ranef(fit2) ), tol=1e-6)
	checkEqualsNumeric( ranef(fit1), unlist(ranef(fit3)), tol=1e-2)


	checkEqualsNumeric( predict(fit1), predict(fit2), tol=1e-7 )
	checkEqualsNumeric( predict(fit1, type="resp"), 
						predict(fit2, type="resp"), 
						tol=1e-6 )

	checkEqualsNumeric( fitted(fit1), exp(fitted(fit2)), tol=1e-6 )

	# checkEqualsNumeric( fixef(fit1), fixef(fit2) )
	checkEqualsNumeric( fixef(fit1), fixef(fit3), tol=1e-3 )
	checkEqualsNumeric( vcov(fit1), vcov(fit2), tol=1e-6  )
	checkEqualsNumeric( coef(summary(fit1))[,1:3], coef(summary(fit2))[,c(1,2,4)], tol=1e-4 )

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

	checkEqualsNumeric( ranef(fit1), unlist(ranef(fit2)), tol=1e-6)
	checkEqualsNumeric( ranef(fit1), unlist(ranef(fit3)), tol=1e-2)


	checkEqualsNumeric( fitted(fit1), exp(fitted(fit2)), tol=1e-6 )
	# checkEqualsNumeric( fixef(fit1), fixef(fit2) )
	checkEqualsNumeric( fixef(fit1), fixef(fit3), tol=1e-3 )
	checkEqualsNumeric( vcov(fit1), vcov(fit2), tol=1e-6  )
	checkEqualsNumeric( coef(summary(fit1))[,1:3], coef(summary(fit2))[,c(1,2,4)], tol=1e-4 )

	fit1$sigSq_g
	fit1$sigSq_e

	fit1 = fastlmm( y ~ x + (1|Indiv), data)
	fit2 = nlme::lme( y ~ x, random = ~ 1|Indiv, data)
	checkEqualsNumeric( sigma(fit1), sigma(fit2), tol=1e-4) 

	# Negative Binomial
	###################

	fit1 = fastglmm.nb( y ~ (1|Indiv), data)
	fit2 = glmer.nb( y ~ (1|Indiv), data)

	family(fit1)
	family(fit2)

	diag(vcov(fit1))
	diag(vcov(fit2))

	coef(summary(fit1))
	coef(summary(fit2))


	fitted(fit1)[1:3]
	fitted(fit2)[1:3]


	w = seq(nrow(data))
	w = w / sum(w)
	# w[] = 1
	fit1 = fastglmm.nb( y ~ (1|Indiv), data, weights=w)
	fit2 = glmer.nb( y ~ (1|Indiv), data, weights=w)

	family(fit1)
	family(fit2)

	diag(vcov(fit1))
	diag(vcov(fit2))

	coef(summary(fit1))
	coef(summary(fit2))


	fitted(fit1)[1:3]
	fitted(fit2)[1:3]




	# library(MASS)
	# library(fastglmm)
	# library(lme4)
	# library(glmmTMB)


	# set.seed(101)
	# dd <- expand.grid(f1 = factor(1:3),
	#                f2 = LETTERS[1:2], g=1:20, rep=1:15,
	#        KEEP.OUT.ATTRS=FALSE)
	# mu <- 5*(-4 + with(dd, as.integer(f1) + 4*as.numeric(f2)))
	# dd$y <- rnbinom(nrow(dd), mu = mu, size =1)
	# dd$g = factor(dd$g)
	# nrow(dd)

	# fit1 <- fastglmm.nb(y ~ f1*f2 + (1|g), data=dd)
	# fit2 <- glmer.nb(y ~ f1*f2 + (1|g), data=dd)
	# fit3 <- fastglmm.nb.glmer(y ~ f1*f2 + (1|g), data=dd)
	# fit4 <- glmer(y ~ f1*f2 + (1|g), data=dd, family=family(fit2))
	# fit5 <- glmmTMB(y ~ f1*f2 + (1|g), data=dd, family=nbinom2)


	# family(fit1)
	# family(fit2)
	# family(fit3)
	# family(fit4)
	# family(fit5)

	# diag(vcov(fit1))
	# diag(vcov(fit2))
	# diag(vcov(fit3))
	# diag(vcov(fit4))

	# coef(summary(fit1))
	# coef(summary(fit2))
	# coef(summary(fit3))
	# coef(summary(fit4))
	# coef(summary(fit5))$cond



	# fitted(fit1)[1:3]
	# fitted(fit2)[1:3]
	# fitted(fit3)[1:3]
	# fitted(fit4)[1:3]



	# n_reps = 10

	# system.time(replicate(n_reps, fastlmm( y ~ x + (1|Indiv), data)))
	# system.time(replicate(n_reps, fastglmm( y ~ x + (1|Indiv), data, family = poisson())))
	# system.time(replicate(n_reps, fastglmm( y ~ x + (1|Indiv), data, family = poisson(), init="glm")))
	# system.time(replicate(n_reps, fastglmm( y ~ x + (1|Indiv), data, family = poisson(), init="fastglm")))
	# system.time(replicate(n_reps, glmmPQL( y ~ x, random = ~ 1|Indiv, data, family = poisson())))
	# system.time(replicate(n_reps, glmer( y ~ x + (1|Indiv), data, family = poisson())))



	# source("/Users/gabrielhoffman/workspace/repos/fastglmm/R/fastglmm.R")

	# fit1 = fastglmm( y ~ x + (1|Indiv), data, family = poisson())


	# library(profvis)
	# profvis({
	# a = replicate(n_reps, fastglmm( y ~ x + (1|Indiv), data, family = poisson()))
	# })



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


	



	# theta_ml2 = function(y, mu, n, weights){

	# 	if( missing(weights) ){
	# 		weights = rep(1, length(y))
	# 	}
	# 	if( missing(n) ){
	# 		n = sum(weights)
	# 	}

	# 	nb_ll = function(y, mu, n, weights, theta){
	# 		tmp = y + theta
	# 		value = lgamma(tmp) - lgamma(theta) - lgamma(y+1) + y*log(mu) + theta * log(theta) - tmp*log(mu + theta)
	# 		crossprod(value, weights)[1] / n
	# 	}
	# 	res = optimize(function(x) -1*nb_ll(y, mu, n, weights, exp(x)), interval=c(-10, 20))
	# 	exp(res$minimum)
	# }

	# w = seq(nrow(data)) / 10.3

	# fit1 = fastglmm.nb( y ~ x + (1|Indiv), data)

	# theta.ml(data$y, fitted(fit1), n = sum(w), weights=w, limit=20300)
	# theta_ml2(data$y, fitted(fit1), n = sum(w), weights=w)



	# system.time(theta.ml(data$y, fitted(fit1), n = sum(w), weights=w, limit=200))
	# system.time(fastglmm:::.theta_ml(data$y, fitted(fit1), n = sum(w), weights=w, limit=200))
	# system.time(theta_ml2(data$y, fitted(fit1)))


	# fastglmm:::.theta_ml(data$y, fitted(fit1), n = sum(w), weights=w, limit=200)

	# fastglmm:::.theta_ml2(data$y, fitted(fit1), n = sum(w), weights=w)




}








