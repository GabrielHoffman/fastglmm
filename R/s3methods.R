
#' ANOVA Tables
#' 
#' ANOVA Tables
#' 
#' @param object fitted model of class \code{fastlmm}
#' @param ... other args, not used
#' 
#' @importFrom dplyr bind_rows mutate `%>%`
#' @importFrom tibble column_to_rownames
#' @export
anova.fastlmm <- function(object,...){

  # Assign each coef to a contrast
  asgn <- attr(object$design, "assign")

  # Names of effects
  nmeffects <- attr(terms(object), "term.labels")[unique(asgn)]

  if( attr(terms(object),"intercept") == 1){
    nmeffects <- c("(Intercept)", nmeffects)
  }

  # Numerator degrees of freedom
  df <- lengths(split(asgn, asgn))

  df1 <- df2 <- NULL

  res = lapply(seq(0, max(asgn)), function(i){

    # create contrast matrix with 1's for this component
    L = rep(0, length(asgn))
    L[which(asgn == i)] = 1

    F.stat <- (L %*% coef(object)) %*% solve(L %*% vcov(object) %*% L) %*% (L %*% coef(object))

    # for comparison
    # res = linearHypothesis(fit, L, test="F")

    data.frame(id = nmeffects[i+1], 
      df1 = sum(L), 
      df2 = df.residual(object),
      F = F.stat)
    }) %>%
    bind_rows %>%
    column_to_rownames("id") %>%
    mutate('Pr(>F)' = pf(F, df1, df2, lower.tail=FALSE))

  structure(res, heading = "Analysis of Variance Table",
                  class = c("anova", "data.frame"))
}

#' ANOVA Tables
#' 
#' ANOVA Tables
#' 
#' @param object fitted model of class \code{fastlmm}
#' @param ... other args, not used
#' 
#' @importFrom dplyr bind_rows mutate `%>%`
#' @importFrom tibble column_to_rownames
#' @export
anova.fastglmm <- function(object,...){

  # Assign each coef to a contrast
  asgn <- attr(object$design, "assign")

  # Names of effects
  nmeffects <- attr(terms(object), "term.labels")[unique(asgn)]

  if( attr(terms(object),"intercept") == 1){
    nmeffects <- c("(Intercept)", nmeffects)
  }

  # Numerator degrees of freedom
  df <- lengths(split(asgn, asgn))

  df1 <- df2 <- NULL

  res = lapply(seq(0, max(asgn)), function(i){

    # create contrast matrix with 1's for this component
    L = rep(0, length(asgn))
    L[which(asgn == i)] = 1

    stat <- (L %*% coef(object)) %*% solve(L %*% vcov(object) %*% L) %*% (L %*% coef(object))

    # for comparison
    data.frame(id = nmeffects[i+1], 
      df = sum(L), 
      Chisq = stat)
    }) %>%
    bind_rows %>%
    column_to_rownames("id") %>%
    mutate('Pr(>Chisq)' = pchisq(Chisq, df, lower.tail=FALSE))

  structure(res, heading = "Analysis of Variance Table",
                  class = c("anova", "data.frame"))
}




#' Extract Model Coefficients
#'
#' Extract Model Coefficients
#'
#' @param object fitted model of class \code{fastlmm}
#' @param ... other args, not used
#'
#' @export
coef.fastlmm <- function(object, ...) {
  object$coefficients
}


#' Cook's Distance Metric
#'
#' Cook's Distance Metric
#'
#' @param model fitted model of class \code{fastlmm}
#' @param ... other args, not used
#'
#' @export
cooks.distance.fastlmm <- function(model,...){

  p <- ncol(model$design)
  h <- hatvalues(model)
  disp <- sigma(model)^2
  resid.pearson <- residuals(model, type="pearson")
  (resid.pearson/(1 - h))^2 * h/(disp * p)
}


#' Model Deviance
#' 
#' Model Deviance
#'
#' @param object fitted model of class \code{fastlmm}
#' @param ... other args, not used
#' 
#' @export
deviance.fastlmm <- function(object, ...) {
  -2*as.numeric(logLik(object))
}

#' @importFrom stats extractAIC
#' @export
extractAIC.fastlmm <- function(fit, scale = 0, k = 2, ...) {
    L <- logLik(fit)
    edf_value <- edf.fastlmm(fit)
    c(edf_value,-2*L + k*edf_value)
}


#' Extract Model Family
#' 
#' Extract Model Family
#'
#' @param object fitted model of class \code{fastlmm}
#' @param ... other args, not used
#' 
#' @importFrom stats family gaussian
#' @rdname family
#' @export
family.fastlmm <- function(object, ...) {
  gaussian()
}


#' @rdname family
#' @export
family.fastglmm <- function(object,...){
    object$family
}


#' Return Diagonals of Hat Matrix
#' 
#' Return Diagonals of Hat Matrix
#'
#' @param model fitted model of class \code{fastlmm}
#' @param ... other args, not used
#' 
#' @importFrom Matrix rowSums
#' @export
hatvalues.fastlmm <- function(model, ...){

  n <- nrow(model$design)
  k <- length(model$s)
  X <- with(model, sqrt(c(weights)) * design)
  Usq <- model$U^2

  h1 <- model$delta*with(model, Usq %*% (1/(s+delta)) + (1 - rowSums(Usq)) / delta)

  A <- with(model, X / delta - U %*% ((s/(delta*s + delta^2)) * crossprod(U, X)))
  D <- solve(crossprod(A, X))
  h2 <- model$delta * rowSums(A * (A %*% D))

  # hatvalues
  1 - h1 + h2
}


#' Effective Degrees-of-Freedom of Model Fit
#' 
#' Effective Degrees-of-Freedom of Model Fit
#'
#' @param object fitted model 
#' @param ... other args, not used
#' 
#' @export
edf <- function(object, ...) {
  UseMethod("edf")
}


#' @export
edf.fastlmm <- function(object, ...){

  n <- nrow(object$design)
  k <- length(object$s)
  X <- with(object, sqrt(c(weights)) * design)

  # sum(h1)
  h1.sum <- object$delta*with(object, sum(1/(s+delta)) + (n-k) / delta)

  # sum(h2)
  A <- with(object, X / delta - U %*% ((s/(delta*s + delta^2)) * crossprod(U, X)))
  D <- solve(crossprod(A, X))
  h2.sum <- object$delta * sum(A * (A %*% D))

  # sum(hatvalues(object))
  n - h1.sum + h2.sum
}


#' Residual Degrees-of-Freedom
#' 
#' Residual Degrees-of-Freedom
#'
#' @param object fitted model of class \code{fastlmm}
#' @param ... other args, not used
#' 
#' @export
df.residual.fastlmm <- function(object, ...){

  n <- nrow(object$design)
  k <- length(object$s)
  X <- with(object, sqrt(c(weights)) * design)

  # sum(h1)
  h1.sum <- object$delta*with(object, sum(1/(s+delta)) + (n-k) / delta)

  # sum(h2)
  A <- with(object, X / delta - U %*% ((s/(delta*s + delta^2)) * crossprod(U, X)))
  D <- solve(crossprod(A, X))
  h2.sum <- object$delta * sum(A * (A %*% D))

  # sum(hatvalues(object))
  # n - h1.sum + h2.sum
  h1.sum - h2.sum
}


#' Test Linear Hypothesis
#' 
#' Test Linear Hypothesis
#' 
#' @param model fitted model of class \code{fastlmm}
#' @param ... other args passed to \code{car::linearHypothesis.default()}
#' 
#' @examples
#' library(MASS)
#' 
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'        family = binomial(), data = bacteria)
#' 
#' linearHypothesis(fit, "trtdrug", test="F")
#
#' @importFrom car linearHypothesis
#' @importFrom car linearHypothesis.default
#' @method linearHypothesis fastlmm
#' 
#' @return \code{car::linearHypothesis()}
#' @seealso \code{car::linearHypothesis()}
#' @rdname linearHypothesis
#' @export
linearHypothesis.fastlmm <- function(model, ...){

  linearHypothesis.default(model, ...)
}

#' @rdname linearHypothesis
#' @export
linearHypothesis <- car::linearHypothesis


#' Extract Log-Likelihood
#' 
#' Extract Log-Likelihood
#'
#' @param object fitted model of class \code{fastlmm}
#' @param ... other args, not used
#'
#' @export
logLik.fastlmm <- function(object, ...) {
  res <- object$residuals
  p <- object$rank
  N <- length(res)
  N0 <- N
  val <- object$logLik
  attr(val, "nall") <- N0
  attr(val, "nobs") <- N
  attr(val, "df") <- edf(object)#p + 2
  class(val) <- "logLik"
  val
}


#' Extract the Number of Observations from a Fit
#' 
#' Extract the Number of Observations from a Fit
#'
#' @param object fitted model of class \code{fastlmm}
#' @param ... other args, not used
#'
#' @importFrom stats nobs
#' @export
nobs.fastlmm <- function(object, ...) {
  if (!is.null(w <- weights(object) ) ) {
    sum(w != 0)
  } else {
    NROW(object$residuals)
  }
}


#' Diagnostic Plots for \code{fastlmm} Fits
#' 
#' Diagnostic Plots for \code{fastlmm} Fits
#' 
#' @param x fitted model of class \code{fastlmm}
#' @param form formula to plot
#' @param ... other args, not used
#' 
#' @importFrom stats loess 
#' @importFrom graphics abline lines 
#' @export
plot.fastlmm <- function(x,
         form = resid(x, type = "pearson") ~ fitted(x),...){

  plot(form, col="dodgerblue")
  abline(h=0)

  # loess smooth
  lw = loess(form)
  j <- order(lw$x)
  lines(lw$x[j],lw$fitted[j],col="red",lwd=3)
}


#' Model Predictions
#'
#' Model Predictions
#'
#' @param object fitted model of class \code{fastlmm}
#' @param newdata optionally, a data frame in which to look for variables with which to predict.  If omitted, the fitted  predictors are used.
#' @param type the type of prediction required.  The default is on the scale of the linear predictors; the alternative \code{"response"} is on the scale of the response variable.  Thus for a default  binomial model the default predictions are of log-odds (probabilities on logit scale) and \code{type = "response"} gives the predicted probabilities.
#' @param ... other args, not used
#'
#' @rdname predict
#' @export
predict.fastlmm <- function(object, newdata = NULL, ...) {
  if (!is.null(newdata)) {
    stop("newdata is not currently supported")
  }
  fitted(object)
}


#' @rdname predict
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


#' Print model
#'
#' Print model
#'
#' @param x fitted model of class \code{fastlmm}
#' @param digits minimal number of _significant_ digits,
#' @param ... other args, not used
#'
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


#' Extract the modes of the random effects
#'
#' Extract the modes of the random effects
#'
#' @param object fitted model of class \code{fastlmm}
#' @param ... other args, not used
#'
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

  # get name of random effect
  id.ranef = findbars(formula(object))

  if( length(id.ranef) > 1){
    stop("Only 1 random effect is supported")
  }
  id.ranef = all.vars(id.ranef[[1]])

  # since U^T U is identity if the GRM is full rank
  v <- with(object, sqrt(s)*ru / (s + delta))
  rownames(v) <- colnames(object$U)
  lst = list()
  lst[[id.ranef]] = as.matrix(v)
  lst
}


#' Extract Studentized Residuals
#' 
#' Extract Studentized Residuals
#' 
#' @param model fitted model of class \code{fastlmm}
#' @param ... other args, not used
#'
#' @importFrom stats fitted
#' @export
rstudent.fastlmm <- function(model, ...) {

  res.dev <- residuals(model, type="deviance")
  h <- hatvalues(model)
  pr <- residuals(model, type="pearson")
  r <- sign(res.dev) * sqrt(res.dev^2 + (h * pr^2)/(1 - h))
  r[is.infinite(r)] <- NaN
  r / sigma(model)
}


#' Extract Residual Standard Deviation 
#' 
#' Extract Residual Standard Deviation 
#' 
#' @param object fitted model of class \code{fastlmm}
#' @param ... other args, not used
#' 
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

  form <- as.character(formula(x))
  cat("Formula: ", paste(form[2], form[1], form[3]), sep='')
  cat("\n\nCoefficients:\n")
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

  cat("\nResidual df:", round(x$rdf, digits=1), "\n")

  cat("\nVariance components:")
  cat("\n  sigSq_g:", format(x$sigSq_g, digits = digits))
  cat("\n  sigSq_e:", format(x$sigSq_e, digits = digits))
  # cat("\n  delta: ", format(x$delta, digits=digits))
  cat("\n  hSq:    ", format(100 * x$sigSq_g / (x$sigSq_g + x$sigSq_e), digits = digits), "%\n")

  cat("\n")
}


#' Object Summaries
#'
#' Object Summaries
#' 
#' @param object fitted model of class \code{fastlmm}
#' @param ... other args, not used
#'
#' @importFrom stats coef pt
#' @export
summary.fastlmm <- function(object, ...) {
  z <- object

  est <- coef(object)
  se <- object$se
  rdf <- df.residual(object)
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
  ans$rdf <- rdf
  ans$sigSq_g <- object$sigSq_g
  ans$sigSq_e <- object$sigSq_e
  ans$delta <- object$delta
  ans$method <- object$method
  ans$formula <- object$formula

  class(ans) <- "summary.fastlmm"
  ans
}


#' Model Terms
#' 
#' Model Terms
#' 
#' @param x fitted model of class \code{fastlmm}
#' @param ... other args, not used
#'
#' @importFrom stats terms
#' @importFrom lme4 nobars
#' @export
terms.fastlmm <- function(x, ...) {
  terms(nobars(formula(x)))
}


#' Calculate Variance-Covariance Matrix for a Fitted Model Object
#' 
#' Calculate Variance-Covariance Matrix for a Fitted Model Object
#' 
#' @param object fitted model of class \code{fastlmm}
#' @param ... other args, not used
#' 
#' @export
vcov.fastlmm <- function(object, ...) {
  object$vcov
}


#' Extract fixed-effects estimates
#' 
#' Extract fixed-effects estimates
#' 
#' @param object fitted model of class \code{fastlmm}
#' @param ... other args, not used
#' 
#' @importFrom lme4 fixef
#' @export
fixef.fastlmm <- function(object, ...) {
  coef(object)
}


#' Extract Model Fitted Values
#' 
#' Extract Model Fitted Values
#' 
#' @param object fitted model of class \code{fastlmm}
#' @param ... other args, not used
#'
#' @rdname fitted
#' @importFrom stats fitted
#' @export
fitted.fastlmm <- function(object, ...) {
  # v <- object$Z %*% ranef.fastlmm(object) + object$design %*% coef(object)

  a <- object$U %*% (sqrt(object$s) * ranef(object)[[1]])
  v <- a / sqrt(object$weights) + object$design %*% coef(object)
  v <- as.numeric(v)

  if (!is.null(object$offset)) {
    v <- v + object$offset
  }
  v
}


# get eta from fastglmm
get_eta = function(object){
  class(object) <- "fastlmm"
  fitted(object)
}

#' @rdname fitted
#' @importFrom stats fitted
#' @export
fitted.fastglmm = function(object,...){

  # convert eta to mu
  object$family$linkinv( get_eta(object) )
}


#' Extract Model Residuals
#' 
#' Extract Model Residuals
#' 
#' @param object fitted model of class \code{fastlmm}
#' @param type the type of residuals which should be returned.  
#' @param ... other args, not used
#'
#' @seealso \code{residuals.glm()}
#' @importFrom stats residuals
#' @export
residuals.fastlmm <- function(object, type = c("working", "response", "deviance", "pearson"), ...) {

  type <- match.arg(type)

  w <- weights(object)
  if( is.null(w) ) w <- 1

  r <- object$y / sqrt(w) - fitted(object)

  if (!is.null(object$offset)) {
    r <- r + object$offset
  }
  r[is.infinite(r)] <- NaN
  
  switch(type, 
    working = , 
    response = r, 
    deviance = , 
    pearson = r * sqrt(w))
}

#' @export
residuals.fastglmm <- function(object, type = c("response", "deviance", "pearson"), ...) {

  type <- match.arg(type)

  if( type %in% c("working", "partial") ){
    stop("Working and partial residuals not supported yet")
  }

  # use prior weights, not PQL weights
  # so set to 1 now
  wts <- object$prior.weights
  if( is.null(wts) ) wts <- 1

  mu <- fitted(object)
  y <- object$response  
  fam <- object$family

  switch(type, 
    deviance = 
      if ( df.residual(object) > 0) {
        d.res <- sqrt(pmax((fam$dev.resids)(y, mu, wts), 0))
        ifelse(y > mu, d.res, -d.res)
      } else rep.int(0, length(mu)), 
    pearson = 
      (y - mu) * sqrt(wts)/sqrt(fam$variance(mu)), 
    # working = r, 
    response = y - mu
    # partial = r
    )
}










