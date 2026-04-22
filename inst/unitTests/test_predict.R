



test_predict = function(){


  library(MASS)
  library(lme4)
  library(fastglmm)
  library(RUnit)

  i = c(3, 45, 219)
  # fit1 <- glmer(y ~ trt + I(week > 2) + (1 | ID),
  #    family = binomial(), data = bacteria[-i,])
  # predict(fit1, newdata=bacteria[i,])

  # IDs already in model
  ######################
  fit2 <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
     family = binomial(), data = bacteria[-i,])
  predict(fit2, newdata=bacteria[i,])

  for( type in c("link", "response")){
    # using new data
    a = predict(fit2, newdata=bacteria, type=type)[i,]
    b = predict(fit2, newdata=bacteria[i,], type=type)
    checkEqualsNumeric(a,b)

    # omitting new data
    a = predict(fit2, type=type)
    b = predict(fit2, newdata=bacteria[-i,], type=type)
    checkEqualsNumeric(a,b)
  }

  # New IDs
  ##########
  i = bacteria$ID == "X01"
  fit2 <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
     family = binomial(), data = bacteria[-i,])
  predict(fit2, newdata=bacteria[i,])

  for( type in c("link", "response")){
    # using new data
    a = predict(fit2, newdata=bacteria, type=type)[i,]
    b = predict(fit2, newdata=bacteria[i,], type=type)
    checkEqualsNumeric(a,b)

    # omitting new data
    a = predict(fit2, type=type)
    b = predict(fit2, newdata=bacteria[-i,], type=type)
    checkEqualsNumeric(a,b)
  }

  # PsychAD
  data(PsychAD)

  form <- PTPRG ~ offset(log(libSize)) + Dx + Age + Sex + (1|SubID)

  keep = PsychAD$SubID != "M10031"
  # fit1 <- glmer(form, PsychAD[keep,], family=poisson(), nAGQ=0)
  fit2 <- fastglmm(form, PsychAD[keep,], family=poisson())

  # a = predict(fit1, newdata=PsychAD[keep,],)
  # b = predict(fit2, newdata=PsychAD[keep,])
  # plot(a,b); abline(0,1, col="red")

  for( type in c("link", "response")){
    # using new data
    a = predict(fit2, newdata=PsychAD, type=type)[which(keep),]
    b = predict(fit2, newdata=PsychAD[keep,], type=type)
    checkEqualsNumeric(a,b)

    # omitting new data
    a = predict(fit2, type=type)
    b = predict(fit2, newdata=PsychAD[keep,], type=type)
    checkEqualsNumeric(a,b)
  }

}