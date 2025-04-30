
test_heritability = function(){


	# IN PROGRESS

	# q()
	# R
	library(fastlmm)
	library(lme4)
	sleepstudy = rbind(sleepstudy, sleepstudy)
	sleepstudy = rbind(sleepstudy, sleepstudy)
	sleepstudy = rbind(sleepstudy, sleepstudy)
	sleepstudy = rbind(sleepstudy, sleepstudy)
	sleepstudy$Subject = sample(LETTERS[1:10], nrow(sleepstudy), replace=TRUE)
	sleepstudy$Subject = factor(sleepstudy$Subject)

	fit <- fastlmm(Reaction ~ Days +  (1 | Subject), sleepstudy)

	fastlmm:::heritability(fit, "info")

	fastlmm:::heritability(fit, "perm", 1000)




}