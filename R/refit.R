
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
#' summary(refit(fit, delta=1000))
#
#' @export
#' @keywords internal
refit = function(fit, delta){

	fastlmm.fit(fit$y, 
		X = fit$design, 
		Z = fit$Z, 
		weights = fit$weights, 
		delta = delta)
}