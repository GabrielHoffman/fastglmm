ll_R <- function(delta, Y, X, Yu, Xu, U, s) {
  info <- QXX <- sig_g <- NA

  n <- nrow(X)
  rank <- nrow(Xu)

  cp_X_low <- crossprod(X) - crossprod(Xu)
  cp_X_low_Y_low <- crossprod(X, Y) - crossprod(Xu, Yu)

  # Eval Beta
  inv_s_delta <- 1 / (s + delta)
  inv_s_delta_Yu <- inv_s_delta * Yu
  inv_s_delta_Xu <- inv_s_delta * Xu

  QXX <- crossprod(Xu, inv_s_delta_Xu) + cp_X_low / delta
  QXY <- crossprod(Xu, inv_s_delta_Yu) + cp_X_low_Y_low / delta
  beta <- solve(QXX, as.matrix(QXY))

  # Eval sig_g
  ru <- Yu - Xu %*% beta
  r <- Y - X %*% beta
  inv_s_delta_ru <- inv_s_delta * ru

  Qrr <- crossprod(ru, inv_s_delta_ru) + (crossprod(r)[1] - crossprod(ru)[1]) / delta
  sig_g <- Qrr[1] / n


  -n / 2 * log(2 * pi * sig_g) - 1 / 2 * (sum(log(s + delta)) + (n - rank) * log(delta)) - n / 2
}

#' Fit linear mixed model using SVD of covariance
#'
#' Fit linear mixed model using SVD of covariance to scale to large sample sizes.
#'
#' @param Y response vector
#' @param X matrix of covariates
#' @param U principal components of covariance matrix
#' @param s eigen values from of covariance matrix
#' @param weights vector weights with value for each sample
#' @param Xu pre-transformed X value
#' @param Yu pre-transformed Y value
#' @param delta ratio of variance components estimated using
#' @param sig_a_fixed if \code{FALSE}, estimate \code{sigSq_a} from data
#' @param rank number of of principal components used
#'
#' @details Fit a linear mixed model with a single variance component.
#'
#' @return summary statistics for model fit, and hypothesis testing
#'
#' @importFrom stats optimize pnorm sd pbeta
#' @export
fastlmm_R <- function(Y, X, U, s, weights = rep(1, nrow(X)), Xu = NULL, Yu = NULL, delta = NULL, sig_a_fixed = FALSE, rank = ncol(U)) {
  rank <- min(rank, ncol(U))

  if (rank < ncol(U)) {
    U <- U[, seq_len(rank), drop = FALSE]
    s <- abs(s[seq_len(rank), drop = FALSE])
  }
  if (is.integer(Y)) {
    Y <- as.numeric(Y)
  }

  if (is.null(Xu)) {
    Xu <- crossprod(U, X)
  }
  if (is.null(Yu)) {
    Yu <- crossprod(U, Y)
  }

  log_interval <- c(10, -10)

  n <- nrow(Y)

  if (is.null(n)) {
    n <- length(Y)
  }

  Gamma_XX <- crossprod(X) - crossprod(Xu)
  Gamma_XY <- crossprod(X, Y) - crossprod(Xu, Yu)

  beta <- sigSq_g <- QXX <- ru <- r <- 1

  i <- 0
  ll <- function(delta_log) {
    i <<- i + 1
    delta <- exp(delta_log)

    # Eval Beta
    inv_s_delta <- 1 / (s + delta)
    inv_s_delta_Yu <- inv_s_delta * Yu
    inv_s_delta_Xu <- inv_s_delta * Xu

    QXX <<- crossprod(Xu, inv_s_delta_Xu) + Gamma_XX / delta
    QXY <- crossprod(Xu, inv_s_delta_Yu) + Gamma_XY / delta
    beta <<- solve(QXX, as.matrix(QXY))

    # Eval sig_g
    ru <<- Yu - Xu %*% beta
    r <<- Y - X %*% beta

    inv_s_delta_ru <- inv_s_delta * ru

    if (sig_a_fixed) {
      sigSq_g <<- 1
    } else {
      Qrr <- crossprod(ru, inv_s_delta_ru) + (crossprod(r)[1] - crossprod(ru)[1]) / delta
      sigSq_g <<- Qrr[1] / n
    }

    -n / 2 * log(2 * pi * sigSq_g) - 1 / 2 * (sum(log(s + delta)) + (n - rank) * log(delta)) - n / 2 + 1 / 2 * sum(log(weights))
  }

  if (is.null(delta)) {
    result <- optimize(ll, log_interval, maximum = TRUE)

    # delta <- result$maximum
    delta <- exp(result$maximum)
  }

  # Need to evaluate ll(), so that obj values are evaluated
  log_L <- ll(delta_log = log(delta))

  beta <- array(beta, dimnames = list(rownames(beta)))
  sigSq_e <- delta * sigSq_g

  ###################################
  # Hypothesis test using Wald test #
  ###################################

  S <- solve(QXX) * sigSq_g
  beta_se <- sqrt(diag(S))

  pValues <- pnorm(abs(beta), 0, beta_se, lower.tail = FALSE) * 2

  df <- sum(s[seq_len(rank)] / (s[seq_len(rank)] + delta))

  res <- list(
    logLik = log_L,
    coefficients = beta,
    se = beta_se,
    vcov = S,
    delta = delta,
    sigSq_g = sigSq_g,
    sigSq_e = sigSq_e,
    iter = i,
    df = df,
    r = r,
    ru = ru,
    pValues = pValues
  )
  class(res) <- "fastlmm"
  return(res)
}
