


test_generics = function(){

  library(tidyverse)
  library(RUnit)
  library(lme4)
  library(fastglmm)

  implemented = c(methods(class = "fastlmm"), 
                  methods(class = "fastglmm")) %>%
                  gsub("\\.fast(g?)lmm", "", .) %>%
                  unique %>%
                  c("formula", "weights") %>%
                  sort

  target = methods(class = "merMod") %>%
                  gsub("\\.merMod", "", .) %>%
                  unique %>%
                  sort

  setdiff(target, implemented)


  f_check_generics = function(fit1, fit2, exclude = c()){

    # check df.residual
    checkEqualsNumeric( n - df.residual(fit1), sum(hatvalues(fit1)))

    tol = 1e-5
    # For each generic function
    exclude = c(exclude, "coef", "print", "plot", "df.residual", "extractAIC", "edf")
    for(fx in setdiff(implemented, exclude) ){

      cat(fx, "\n")

      if( fx == "linearHypothesis"){

        # run on fastlmm fit
        res1 = get(fx)( fit1, "Days")

        # run on lmer fit
        res2 = get(fx)( fit2, "Days")

      }else{
        # run on fastlmm fit
        res1 = get(fx)( fit1 )

        # run on lmer fit
        res2 = get(fx)( fit2 )
      }

      if( fx %in% c("family")){
        checkEquals( res1, res2 )
      }else if( fx %in% c("formula")){
        identical(res1, res2)
      }else if( fx %in% c("anova")){
        checkEqualsNumeric(anova(fit1)$F[-1], anova(fit2)$`F value`, tol=1e-6)
      }else if( fx %in% c("linearHypothesis")){
        checkEqualsNumeric(res1$Chisq, res2$Chisq, tol=1e-6)
      }else if( fx %in% c("terms")){
        ids = intersect(names(attributes(res1)), names(attributes(res2)))
        checkEquals(attributes(res1)[ids], attributes(res2)[ids])
      }else if( fx %in% c("ranef")){
        checkEqualsNumeric( res1[[1]], res2[[1]], tol=tol)
      }else if( fx %in% c("weights")){
        checkEqualsNumeric( c(res1), c(res2), tol=tol)
      }else if( fx %in% c("residuals")){

        types <- c("deviance", "pearson", "working","response")
        for(type in types){
          res1 = get(fx)( fit1, type)
          res2 = get(fx)( fit2, type)

          checkEqualsNumeric( res1, res2, tol=tol )
        }

      }else if( fx %in% c("summary")){
        checkEqualsNumeric( coef(res1)[,1:3], coef(res2)[,1:3], tol=tol)
      }else{
        checkEqualsNumeric( res1, res2, tol=tol )
      }
    }
  }



  # Weighted LMM
  w = seq(nrow(sleepstudy))
  w = w / mean(w)
  n = nrow(sleepstudy)
  form = Reaction ~ Days + (1 | Subject)

  fit1 <- fastlmm(form, sleepstudy, weights = w)
  fit2 <- lmer(form, sleepstudy, REML = FALSE, weights = w)
  f_check_generics(fit1, fit2)

  # Weighted LMM with offset
  set.seed(1)
  sleepstudy$off = rpois(nrow(sleepstudy), 1000)
  form = Reaction ~ Days + (1 | Subject) + offset(log(off))

  fit1 <- fastlmm(form, sleepstudy, weights = w)
  fit2 <- lmer(form, sleepstudy, REML = FALSE, weights = w)
  f_check_generics(fit1, fit2)

  fit3 <- fastglmm(form, sleepstudy, weights = w)
  f_check_generics(fit1, fit3, exclude = "anova")




}
























