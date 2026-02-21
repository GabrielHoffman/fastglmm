

# https://chatgpt.com/share/698f66b1-88b0-800b-9001-e28f444f6015
est_hessian <- function(fit){

  delta <- fit$delta
  sigSq_g <- fit$sigSq_g
  s <- fit$s

  n <- length(fit$y)
  r <- length(s)
  H <- matrix(0, 2,2)

  a <- sum(1/(s + delta)) + (n-r)/delta
  b <- sum(1/(s + delta)^2) + (n-r)/delta^2

  H[1,1] <- (n - 2*delta*a + delta^2*b) / (2*sigSq_g^2)
  H[1,2] <- H[2,1] <- (a - delta*b) / (2*sigSq_g^2)
  H[2,2] <- b / (2*sigSq_g^2)
  H
}



est_gradient <- function(fit, L){

  # slow versions
  # W <- with(fit, solve(tcrossprod(Z) + diag(delta, nrow(Z))))
  # A <- crossprod(X, W) %*% X
  # B <- crossprod(X, W %*% W) %*% X

  Xu <- crossprod(fit$U, fit$design)
  Gamma_XX <- crossprod(fit$design) - crossprod(Xu)
  inv_s_delta <- 1 / (fit$s + fit$delta)
  inv_s_delta_Xu <- inv_s_delta * Xu

  A <- crossprod(Xu, inv_s_delta_Xu) + Gamma_XX / fit$delta

  inv_s_delta_Xu <- inv_s_delta^2 * Xu
  B <- crossprod(Xu, inv_s_delta_Xu) + Gamma_XX / fit$delta^2

  lapply(seq(nrow(L)), function(i){

    g <- c(0, 0)
    invAL <- solve(A, L[i,])
    C <- invAL %*% B %*% invAL
    g[1] <- crossprod(L[i,], invAL) - fit$delta*C
    g[2] <- C
    g
  })
}

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

  H <- est_hessian(fit)
  g <- est_gradient(fit, L)

  sapply(seq(nrow(L)), function(i){
    var_Lbeta <- crossprod(L[i,], vcov(fit)) %*% L[i,]
    v_numerator <- 2 * var_Lbeta^2
    v_denom <- crossprod(g[[i]], solve(H)) %*% g[[i]]

    v_numerator / v_denom  
  })
}

