# # 1) create profile log-likelihood surface from delta surface, and get standard error from hessian of hsq surface at hsq_hat
# # 2) fast permutations

# # emdbook::pchibarsq
# #' @importFrom stats pchisq
# pchibarsq <- function (p, df = 1, mix = 0.5, lower.tail = TRUE, log.p = FALSE) {
#     df <- rep(df, length.out = length(p))
#     mix <- rep(mix, length.out = length(p))
#     c1 <- ifelse(df == 1, if (lower.tail) 
#         1
#     else 0, pchisq(p, df - 1, lower.tail = lower.tail))
#     c2 <- pchisq(p, df, lower.tail = lower.tail)
#     r <- mix * c1 + (1 - mix) * c2
#     if (log.p) 
#         log(r)
#     else r
# }


# #' Hypothesis test for heritability
# #' 
# #' Hypothesis test for heritability using either asymptotic normal approximation  or empirical permutations.
# #' 
# #' @param fit model fit of class \code{fastlmm}
# #' @param method \code{"information"} or \code{"permutation"}, 
# #' @param nperms number of permutations
# #' 
# #' @details For \code{method == "information"}, the profile log-likelihood is evaluted with respect to hsq.  The standard error is obtained from the information matrix based on the Hessian evaluated at the MLE of hsq. The p-value is then computed from this estimated standard error using a normal approximation.  This approach can perform well for large sample sizes, by performs poorly for moderate sample sizes.
# #' 
# #' For \code{method == "permutation"},  ....

# #' @references
# #' Abney, M. (2015). Permutation testing in the presence of polygenic variation. Genetic epidemiology, 39(4), 249-258. \doi{10.1002/gepi.21893}
# #' 
# #' @importFrom numDeriv hessian
# #' @export
# randomEffect.test <- function(fit, method = c("information", "permutation"), nperms = 100) {

#   method <- match.arg(method)

#   # estimate of hsq given delta
#   # hsq_hat <- 1 - 1 / (1 + 1 / fit$delta)
#   hsq_hat <- 1 / (fit$delta + 1)
#   # with(fit, sig_g / (sig_g + sig_e))

#   if (method == "information") {
#     Yu <- crossprod(fit$U, fit$y)
#     Xu <- crossprod(fit$U, fit$design)
#     f <- function(hsq) {
#       # delta <- 1 / (1 / hsq - 1)
#       delta <- hsq / (1 - hsq)
#       -1 * fastglmm:::ll_R(delta, fit$y, fit$design, Yu, Xu, fit$U, fit$s)
#     }

#     # using refit
#     # which is faster for large data?
#     # g <- function(hsq) {
#     #   delta <- 1 / (1 / hsq - 1)
#     #   fit2 <- refit(fit, delta=delta)
#     #   -1*fit2$logLik
#     # }

#     # Fisher information at delta_hat
#     infor <- hessian(f, hsq_hat)

#     # variance of estimate
#     se_hsq <- sqrt(1 / infor)

#     stat <- hsq_hat / se_hsq
#     p.value <- pchisq(stat^2, 1, lower.tail=FALSE)

#     res <- data.frame(hsq = hsq_hat, se = se_hsq, p.value)

#     # based on hessian of sigSq_g
    
  

#   } else if (method == "permutation") {

#     fit_null <- refitModel(fit, delta = 1e-4)
#     residValues <- residuals(fit_null)

#     # Use residuals instead here??
#     fitList <- lapply(seq(nperms), function(i) {
#       y.perm <- sample(fit$y, length(fit$y), replace = TRUE)
#       # r <- sample(residValues, length(residValues), replace = TRUE)
#       # r + predict(fit)

#       fastlmm.fit(y.perm, fit$design, Z = fit$Z)
#     })

#     h_sq_null <- sapply(fitList, function(fit) {
#       # 1 - 1 / (1 + 1 / fit$delta)
#       1 / (fit$delta + 1)
#     })

#     # Approximate null distribution with beta
#     mu <- mean(h_sq_null)
#     se <- sd(h_sq_null)
#     alpha <- mu * (mu * (1 - mu) / se^2 - 1)
#     beta <- (1 - mu) * (mu * (1 - mu) / se^2 - 1)

#     # method of moments for beta distribiton
#     p.value <- pbeta(hsq_hat, alpha, beta, lower.tail = FALSE)

#     res <- data.frame(hsq = hsq_hat, se = se, p.value)
#   }

#   res
# }

