# methods(class="fastlmm")
# getS3method("deviance", "lm")
# getS3method("family", "lm")
# getS3method("print", "summary.lm")

# summary(lm.D9)

# blup

# Extract Model Coefficients
#
# Extract Model Coefficients
#
# @param object fitted model of class \code{fastlmm}
# @param ... other args, not used
#
#' @export
coef.fastlmm <- function(object, ...) {
  object$coefficients
}

# confint(fit2)

# cooks.distance.fastlmm


#' @importFrom stats family gaussian
#' @export
family.fastlmm <- function(object, ...) {
  gaussian()
}

# formula.fastlmm

# hatvalues.fastlmm

# influence.fastlmm

# kappa.fastlmm

# labels.fastlmm

#' @export
logLik.fastlmm <- function(object, ...) {
  res <- object$residuals
  p <- object$rank
  N <- length(res)
  N0 <- N
  val <- object$logLik
  attr(val, "nall") <- N0
  attr(val, "nobs") <- N
  attr(val, "df") <- p + 2
  class(val) <- "logLik"
  val
}

# model.frame.fastlmm

# model.matrix.fastlmm

#' @importFrom stats nobs
#' @export
nobs.fastlmm <- function(object, ...) {
  if (!is.null(w <- object$weights)) {
    sum(w != 0)
  } else {
    NROW(object$residuals)
  }
}

# plot.fastlmm


#' @importFrom stats nobs
#' @export
predict.fastlmm <- function(object, newdata = NULL, ...) {
  if (!is.null(newdata)) {
    stop("newdata is not currently supported")
  }
  fitted(object)
}


#' @importFrom stats coef
#' @export
print.fastlmm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(attr(x, "call")), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2L, quote = FALSE)
  } else {
    cat("No coefficients\n")
  }
  cat("\n")
  invisible(x)
}

# ranef.fastlmm



# rstandard.fastlmm

# rstudent.fastlmm


# show.fastlmm


#' @importFrom stats sigma
#' @export
sigma.fastlmm <- function(object, ...) {
  sqrt(object$sigSq_e)
}

#' @importFrom stats printCoefmat
#' @export
print.summary.fastlmm <- function(
    x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor,
    signif.stars = getOption("show.signif.stars"), ...) {
  if (x$method %in% c("ML", "REML")) {
    cat("Linear mixed model fit by", x$method, "['fastlmm']\n")
  } else {
    cat("Generlized linear mixed model fit by", x$method, "['fastglmm']\n")
  }

  cat("\nCoefficients:\n")
  coefs <- x$coefficients
  if (any(aliased <- x$aliased)) {
    cn <- names(aliased)
    coefs <- matrix(NA, length(aliased), 4, dimnames = list(
      cn,
      colnames(coefs)
    ))
    coefs[!aliased, ] <- x$coefficients
  }
  printCoefmat(coefs, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)

  cat("\nVariance components:")
  cat("\n  sigSq_g:", format(x$sigSq_g, digits = digits))
  cat("\n  sigSq_e:", format(x$sigSq_e, digits = digits))
  # cat("\n  delta: ", format(x$delta, digits=digits))
  cat("\n  hSq:    ", format(100 * x$sigSq_g / (x$sigSq_g + x$sigSq_e), digits = digits), "%\n")

  cat("\n")
}

#' @importFrom stats coef pt
#' @export
summary.fastlmm <- function(object, ...) {
  z <- object

  est <- coef(object)
  se <- object$se
  rdf <- object$df.residual
  tval <- est / se
  ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
  ans$aliased <- is.na(coef(object))
  ans$residuals <- z$residuals
  ans$coefficients <- cbind(
    Estimate = est, `Std. Error` = se,
    `t value` = tval, `Pr(>|t|)` = 2 * pt(abs(tval), rdf,
      lower.tail = FALSE
    )
  )
  ans$sigSq_g <- object$sigSq_g
  ans$sigSq_e <- object$sigSq_e
  ans$delta <- object$delta
  ans$method <- object$method

  class(ans) <- "summary.fastlmm"
  ans
}


#' @export
vcov.fastlmm <- function(object, ...) {
  object$vcov
}

#' @importFrom lme4 ranef
#' @importFrom Matrix crossprod
#' @export
ranef.fastlmm <- function(object, ...) {

  # original, uses Z
  # Zw <- c(sqrt(object$weights)) * object$Z
  # A <- crossprod(object$U, Zw)
  # b <- object$ru / (object$s + object$delta)
  # v <- crossprod(A, b)
  # as.matrix(v)

  U <- s <- weights <- ru <- delta <- NULL

  # Use U and s, but not Z
  # since Z = U diag(sqrt(s))
  # A <- with(object, 
  #       crossprod(U, U %*% Diagonal(length(s), sqrt(s))))
  # b <- with(object, ru / (s + delta))
  # v <- crossprod(A, b)
  # rownames(v) <- colnames(object$Z)
  # as.matrix(v)

  # since U^T U is identity if the GRM is full rank
  v <- with(object, sqrt(s)*ru / (s + delta))
  rownames(v) <- colnames(object$U)
  as.matrix(v)
}

# attempt to simply fitted()
# object$Z %*% v

# with(object, (U * sqrt(s)) %*% crossprod(crossprod(U, c(sqrt(weights)) * U * sqrt(s)), b)) 

# U = object$U
# s = object$s
# Us = (U * sqrt(s))
# weights = object$weights
# Us %*% crossprod(crossprod(U, c(sqrt(weights)) * Us), b) 

#' @importFrom lme4 fixef
#' @export
fixef.fastlmm <- function(object, ...) {
  coef(object)
}

#' @importFrom stats fitted
#' @export
fitted.fastlmm <- function(object, ...) {
  # v <- object$Z %*% ranef.fastlmm(object) + object$design %*% coef(object)

  a <- object$U %*% (sqrt(object$s) * ranef.fastlmm(object))
  v <- a / sqrt(object$weights) + object$design %*% coef(object)
  v <- as.numeric(v)

  if (!is.null(object$offset)) {
    v <- v + object$offset
  }
  v
}


# See ?residuals.glm
#' @importFrom stats residuals
#' @export
residuals.fastlmm <- function(object, ...) {
  if (is.null(object$weights)) {
    w <- 1
  } else {
    w <- object$weights
  }

  v <- object$y / sqrt(w) - fitted(object)

  if (!is.null(object$offset)) {
    v <- v + object$offset
  }
  v
}


# fastglmm
##########

# get eta from fastglmm
get_eta = function(object){
  class(object) <- "fastlmm"
  fitted(object)
}

#' @importFrom stats fitted
#' @export
fitted.fastglmm = function(object,...){

  # convert eta to mu
  object$family$linkinv( get_eta(object) )
}


#' @importFrom stats nobs
#' @export
predict.fastglmm <- function(object, newdata = NULL,
                 type = c("link", "response"), ...) {

  type <- match.arg(type)

  if (!is.null(newdata)) {
    stop("newdata is not currently supported")
  }
  
  switch(type,
      link = get_eta(object),
      response = fitted(object))
}

