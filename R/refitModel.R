
#' Refit fastlmm model with new delta
#' 
#' Refit fastlmm model with new delta value
#' 
#' @param fit model fit of class \code{fastlmm}
#' @param delta new value for delta 
#' 
#' @details Useful for evaluating log-likelihood and multiple values of delta
#' 
#' @examples
#' library(lme4)
#' 
#' fit <- fastlmm(Reaction ~ Days + (1 | Subject), sleepstudy)
#' 
#' summary(fit)
#' 
#' summary(refitModel(fit, delta=1000))
#
#' @export
refitModel <- function(fit, delta){
  if (is(fit$U, "sparseMatrix")) {
    fxn <- .fastlmm_vms
  } else {
    fxn <- .fastlmm_vmm
  }

  REML <- FALSE

  res <- fxn(
      y = fit$y,
      X = fit$design,
      U = fit$U,
      s = fit$s,
      weights = weights(fit),
      REML = REML,
      delta = delta,
      left = 1, 
      right = 1, 
      tol = 1e-4, 
      nthreads = 1)

  os = 1
  res <- as.fastlmm(res, design = fit$design, offset = os, method = ifelse(REML, "REML", "ML"))
  res$U <- fit$U
  res$s <- fit$s
  res$formula <- formula(fit)

  res
}
