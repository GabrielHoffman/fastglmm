



test_varpart_offset = function(){

  library(fastglmm)
  library(MASS)
  library(RUnit)

  data(PsychAD)

  X = model.matrix(~Age, PsychAD)
  eta = log(PsychAD$libSize) + X %*% c(0,.02)

  PsychAD$PTPRG = rnegbin(nrow(PsychAD), exp(eta), theta=4)

  # NB model
  form <- PTPRG ~ offset(log(libSize)) + Dx + Age + Sex + (1|SubID)
  fit <- fastglmm.nb(form, PsychAD)
  vp1 = varpart(fit)

  # log fractions
  form <- log(PTPRG+1) - log(libSize) ~ Dx + Age + Sex + (1|SubID)
  fit <- fastlmm(form, PsychAD)
  vp2 = varpart(fit)

  plot(vp1[-5], vp2)

  checkIdentical(max(abs(vp1[-5] - vp2)) < 1e-3, TRUE)


}
