

#' Denominator Degrees of Freedom
#'
#' Denominator degrees of freedom using Satterthwaite method
#' 
#' @param fit model fit from \code{fastlmm()} or \code{fastglmm()}
#' @param L matrix of coefficient contrasts, one per _row_
#'
#' @return array, denominator degrees of freedom for each contrast (i.e. _row_)
#' 
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#' 
#' # denominator degrees of freedom 
#' # used for hypothesis testing below
#' ddf(fit)
#' 
#' # summarize fit
#' summary(fit)
#
#' @export 
ddf <- function(fit, L = diag(1, length(coef(fit))) ){

  stopifnot(is(fit, "fastlmm"))

  .ddf(fit$delta, fit$hessian.vc, fit$A.sat, fit$B.sat, vcov(fit), L )
}

#' @keywords internal
#' @export
.ddf <- function(delta, hessian.vc, A.sat, B.sat, V, L){

  sapply(seq(nrow(L)), function(i){

    invAL <- solve(A.sat, L[i,])
    C <- invAL %*% B.sat %*% invAL
    g <- c(crossprod(L[i,], invAL) - delta*C, C)

    var_Lbeta <- crossprod(L[i,], V %*% L[i,])
    v_numerator <- 2 * var_Lbeta^2
    v_denom <- solve(hessian.vc, g) %*% g

    v_numerator / v_denom  
  })
}

#' DDF for joint test
#'
#' DDF for joint test by summarizing ddf for each coefficient
#'
#' @param nu array of ddf values
#' @param tol tolerance
#'
#' @seealso \code{lmerTest:::get_Fstat_ddf()}
#' @keywords internal
#' @export
get_Fstat_ddf <- function(nu, tol=1e-8) {
  # Computes denominator df for an F-statistic that is derived from a sum of
  # squared t-statistics each with nu_m degrees of freedom.
  #
  # nu : vector of denominator df for the t-statistics
  # tol: tolerance on the consequtive differences between elements of nu to
  #      determine if mean(nu) should be returned.
  #
  # Result: a numeric scalar
  #
  # Returns nu if length(nu) == 1. Returns mean(nu) if all(abs(diff(nu)) < tol;
  # otherwise ddf appears to be downward biased.
  fun <- function(nu) {
    if(any(nu <= 2)) 2 else {
      E <- sum(nu / (nu - 2))
      2 * E / (E - (length(nu))) # q = length(nu) : number of t-statistics
    }
  }
  stopifnot(length(nu) >= 1,
            # all(nu > 0), # returns 2 if any(nu < 2)
            all(sapply(nu, is.numeric)))
  if(length(nu) == 1L) return(nu)
  if(all(abs(diff(nu)) < tol)) return(mean(nu))
  if(!is.list(nu)) fun(nu) else vapply(nu, fun, numeric(1L))
}




# # https://chatgpt.com/share/698f66b1-88b0-800b-9001-e28f444f6015
# est_hessian <- function(fit){

#   delta <- fit$delta
#   sigSq_g <- fit$sigSq_g
#   s <- fit$s

#   n <- length(fit$y)
#   r <- length(s)
#   H <- matrix(0, 2,2)

#   a <- sum(1/(s + delta)) + (n-r)/delta
#   b <- sum(1/(s + delta)^2) + (n-r)/delta^2

#   H[1,1] <- (n - 2*delta*a + delta^2*b) / (2*sigSq_g^2)
#   H[1,2] <- H[2,1] <- (a - delta*b) / (2*sigSq_g^2)
#   H[2,2] <- b / (2*sigSq_g^2)
#   H
# }



# est_gradient <- function(fit, L){

#   # slow versions
#   # W <- with(fit, solve(tcrossprod(Z) + diag(delta, nrow(Z))))
#   # A <- crossprod(X, W) %*% X
#   # B <- crossprod(X, W %*% W) %*% X

#   Xu <- crossprod(fit$U, fit$design)
#   Gamma_XX <- crossprod(fit$design) - crossprod(Xu)
#   inv_s_delta <- 1 / (fit$s + fit$delta)
#   inv_s_delta_Xu <- inv_s_delta * Xu

#   A <- crossprod(Xu, inv_s_delta_Xu) + Gamma_XX / fit$delta

#   inv_s_delta_Xu <- inv_s_delta^2 * Xu
#   B <- crossprod(Xu, inv_s_delta_Xu) + Gamma_XX / fit$delta^2

#   lapply(seq(nrow(L)), function(i){

#     g <- c(0, 0)
#     invAL <- solve(A, L[i,])
#     C <- invAL %*% B %*% invAL
#     g[1] <- crossprod(L[i,], invAL) - fit$delta*C
#     g[2] <- C
#     g
#   })
# }


# ddf.orig <- function(fit, L = diag(1, length(coef(fit))) ){

#   stopifnot(is(fit, "fastlmm"))

#   H <- est_hessian(fit)
#   g <- est_gradient(fit, L)

#   sapply(seq(nrow(L)), function(i){
#     var_Lbeta <- crossprod(L[i,], vcov(fit)) %*% L[i,]
#     v_numerator <- 2 * var_Lbeta^2
#     v_denom <- crossprod(g[[i]], solve(H)) %*% g[[i]]

#     v_numerator / v_denom  
#   })
# }
