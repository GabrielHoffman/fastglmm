
#' ANOVA Tables
#' 
#' ANOVA Tables
#' 
#' @param object fitted model of class \code{fastlmm}
#' @param ... other args, not used
#'
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#'
#' anova(fit)
#
#' @rdname anova
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


#' @rdname anova
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

  df1 <- df2 <- Chisq <- NULL

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
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#'
#' coef(fit)
#
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
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#'
#' cooks.distance(fit)[1:3]
#
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
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#'
#' deviance(fit)
#
#' @export
deviance.fastlmm <- function(object, ...) {
  -2*as.numeric(logLik(object))
}

#' Residual Degrees-of-Freedom
#' 
#' Residual Degrees-of-Freedom
#'
#' @param object fitted model of class \code{fastlmm}
#' @param ... other args, not used
#'
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#'
#' df.residual(fit)
#
#' @export
df.residual.fastlmm <- function(object, ...){

  n <- nrow(object$design)
  k <- length(object$s)
  X <- with(object, sqrt(c(weights)) * design)

  # sum(h1)
  h1.sum <- with(object, delta*sum(1/(s+delta))) + (n-k)

  # sum(h2)
  A <- with(object, X / delta - U %*% ((s/(delta*s + delta^2)) * crossprod(U, X)))
  D <- solve(crossprod(A, X))
  h2.sum <- object$delta * sum(A * (A %*% D))

  # sum(hatvalues(object))
  # n - h1.sum + h2.sum
  h1.sum - h2.sum
}

#' Extract AIC from a Fitted Model
#'
#' Extract AIC from a Fitted Model
#'
#' @param fit fitted model of class \code{fastlmm}
#' @param scale not used
#' @param k numeric specifying the \code{weight} of the _equivalent degrees of freedom_ (=: \code{edf}) part in the AIC formula.
#' @param ... other args, not used
#'
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#'
#' extractAIC(fit)
#
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
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#'
#' family(fit)
#
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


#' @rdname family
#' @export
family.glmmPQL <- function(object,...){
    object$family
}


#' Return Diagonals of Hat Matrix
#' 
#' Return Diagonals of Hat Matrix
#'
#' @param model fitted model of class \code{fastlmm}
#' @param ... other args, not used
#' 
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#'
#' hatvalues(fit)[1:3]
#
#' @importFrom Matrix rowSums
#' @export
hatvalues.fastlmm <- function(model, ...){

  n <- nrow(model$design)
  k <- length(model$s)
  X <- with(model, sqrt(c(weights)) * design)
  Usq <- model$U^2

  h1 <- model$delta*with(model, Usq %*% (1/(s+delta))) + (1 - rowSums(Usq))

  A <- with(model, X / delta - U %*% ((s/(delta*s + delta^2)) * crossprod(U, X)))
  D <- solve(crossprod(A, X))
  h2 <- model$delta * rowSums(A * (A %*% D))

  # hatvalues
  1 - as.numeric(h1) + h2
}


#' Effective Degrees-of-Freedom of Model Fit
#' 
#' Effective Degrees-of-Freedom of Model Fit
#'
#' @param object fitted model 
#' @param ... other args, not used
#' 
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#'
#' edf(fit)
#
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






#' Extract Model Fitted Values
#' 
#' Extract Model Fitted Values
#' 
#' @param object fitted model of class \code{fastlmm}
#' @param ... other args, not used
#'
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#'
#' fitted(fit)[1:3]
#
#' @rdname fitted
#' @importFrom stats fitted
#' @importFrom Matrix crossprod
#' @export
fitted.fastlmm <- function(object, ...) {
    
  blp = ranef(object)[[1]]
  a <- with(object, U %*% (sqrt(s) * crossprod(V, blp)))
  v <- as.matrix(a) / sqrt(object$weights) + object$design %*% coef(object)

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
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
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
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#'
#' logLik(fit)
#
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

#' Family function for Negative Binomial GLMs
#' 
#' Specifies the information required to fit a Negative Binomial generalized linear model, with known \code{theta} parameter.  
#'
#' @param theta The known value of the additional parameter, \code{theta}
#' @param link The link function, as a character string, name or one-element character vector specifying one of \code{log}, \code{sqrt} or \code{identity}, or an object of class \code{"link-glm"}.
#'
#' @export 
negative.binomial <- MASS::negative.binomial

#' Extract the Number of Observations from a Fit
#' 
#' Extract the Number of Observations from a Fit
#'
#' @param object fitted model of class \code{fastlmm}
#' @param ... other args, not used
#'
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#'
#' nobs(fit)
#
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
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#'
#' plot(fit)
#
# @importFrom stats loess 
# @importFrom graphics abline lines 
#' @importFrom graphics abline 
#' @export
plot.fastlmm <- function(x,
         form = resid(x, type = "pearson") ~ fitted(x),...){

  plot(form, col="dodgerblue")
  abline(h=0)

  # too slow for large datasets
  # loess smooth
  # lw = loess(form)
  # j <- order(lw$x)
  # lines(lw$x[j],lw$fitted[j],col="red",lwd=3)
}


#' Model Predictions
#'
#' Model Predictions
#'
#' @param object fitted model of class \code{fastlmm}
#' @param newdata optionally, a data frame in which to look for variables with which to predict.  If omitted, the fitted  predictors are used.
#' @param type the type of prediction required.  The default is on the scale of the linear predictors; the alternative \code{"response"} is on the scale of the response variable.  Thus for a default  binomial model the default predictions are of log-odds (probabilities on logit scale) and \code{type = "response"} gives the predicted probabilities. \code{type = "terms"} computes the linear predict for each model term
#' @param ... other args, not used
#'
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#'
#' predict(fit)[1:3]
#
#' @rdname predict
#' @export
predict.fastlmm <- function(object, newdata = NULL, type = c("response", "terms"), ...) {

  type <- match.arg(type)

  if (!is.null(newdata)) {
    stop("newdata is not currently supported")
  }

  switch(type, 
    response = fitted(object),
    terms = predict_terms(object))
}


#' @rdname predict
#' @export
predict.fastglmm <- function(object, newdata = NULL,
                 type = c("link", "response", "terms"), ...) {

  type <- match.arg(type)

  if (!is.null(newdata)) {
    stop("newdata is not currently supported")
  }
  
  switch(type,
      link = get_eta(object),
      response = fitted(object),
      terms = predict_terms(object))
}


# Compute eta value for each fixed effect
predict_terms <- function( object ){

  # Assign each coef to a contrast
  asgn <- attr(object$design, "assign")

  # Names of effects
  nmeffects <- attr(terms(object), "term.labels")[unique(asgn)]

  if( attr(terms(object),"intercept") == 1){
    nmeffects <- c("(Intercept)", nmeffects)
  }

  # prediction for each term
  Eta <- lapply(seq(0, max(asgn)), function(i){
    # indeces in this term
    idx <- which(asgn == i)
    eta_i <- object$design[,idx,drop=FALSE] %*% coef(object)[idx,drop=FALSE]
    c(eta_i)
    }) 
  Eta <- do.call("cbind", Eta)
  colnames(Eta) <- nmeffects

  Eta
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
#' @keywords internal
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




#' Extract Studentized Residuals
#' 
#' Extract Studentized Residuals
#' 
#' @param model fitted model of class \code{fastlmm}
#' @param ... other args, not used
#'
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#'
#' rstudent(fit)[1:3]
#
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



#' Extract the modes of the random effects
#'
#' Extract the modes of the random effects
#'
#' @param object fitted model of class \code{fastlmm}
#' @param ... other args, not used
#'
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#'
#' ranef(fit)
#
#' @importFrom lme4 ranef
#' @importFrom Matrix crossprod
#' @method ranef fastlmm
#'
#' @rdname ranef
#' @export
ranef.fastlmm <- function(object, ...) {

  # get name of random effect
  id.ranef <- findbars(formula(object))

  if( length(id.ranef) > 1){
    stop("Only 1 random effect is supported")
  }
  id.ranef <- all.vars(id.ranef[[1]])

  weights <- U <- s <- ru <- delta <- NULL

  # Uses Z
  Zw <- with( object, sqrt(c(weights)) * Z)
  A <- crossprod(object$U, Zw)
  b <- with(object, ru / (s + delta))
  v <- as.matrix(crossprod(A, b))
  rownames(v) <- colnames(object$Z)

  # Use U and s, but not Z
  # since Z = U diag(sqrt(s))
  # A <- with(object, 
  #       crossprod(U, U %*% Diagonal(length(s), sqrt(s))))
  # b <- with(object, ru / (s + delta))
  # v <- crossprod(A, b)
  # rownames(v) <- colnames(object$Z)
  # as.matrix(v)
  
  # NEEDS to use V for ordering columns
  # since U^T U is identity if the GRM is full rank
  # v <- with(object, sqrt(s)*ru / (s + delta))

  lst = list()
  lst[[id.ranef]] <- v
  lst
}



#' Extract Fixed Effects
#'
#' Extract Fixed Effects
#'
#' @param object fitted model of class \code{fastlmm}
#' @param ... other args, not used
#'
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#'
#' fixef(fit)
#
#' @importFrom lme4 fixef
#' @method fixef fastlmm
#'
#' @export
fixef.fastlmm <- function(object,...){
  coef(object)
}


#' Extract Residual Standard Deviation 
#' 
#' Extract Residual Standard Deviation 
#' 
#' @param object fitted model of class \code{fastlmm}
#' @param ... other args, not used
#' 
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#'
#' sigma(fit)
#
#' @importFrom stats sigma
#' @export
sigma.fastlmm <- function(object, ...) {
  sqrt(object$sigSq_e)
}


cat.f <- function(...) cat(..., fill = TRUE)

.prt.family <- function(famL) {
  if (!is.null(f <- famL$family)) {
    cat.f(" Family:", f,
        if(!is.null(ll <- famL$link)) paste(" (", ll, ")"))
  }
}

#' @importFrom stats printCoefmat family
#' @export
print.summary.fastlmm <- function(
    x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor,
    signif.stars = getOption("show.signif.stars"), ...) {
  if (x$method %in% c("ML", "REML")) {
    cat("Linear mixed model fit by", x$method, "['fastlmm']\n")
  } else {
    cat("Generlized linear mixed model fit by", x$method, "['fastglmm']\n")
    .prt.family( x$family )
  }

  form <- as.character(formula(x))
  cat(" Formula: ", paste(form[2], form[1], form[3]), sep='')
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

  # cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ", format(x$dispersion), ")\n")

  cat("\nResidual df:", round(x$rdf, digits=1), "\n")

  cat("\nVariance components:")
  cat("\n  sigSq_g:", format(x$sigSq_g, digits = digits))
  cat("\n  sigSq_e:", format(x$sigSq_e, digits = digits))
  cat("\n  delta:  ", format(x$delta, digits=digits))

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

  # if object a GLMM and the family is _not_ gaussian
  isGaussian <- (family(object)$family == "gaussian")
  if( is(object, "fastglmm") & !isGaussian){
    # Normal approximation
    ans$coefficients <- cbind(
      Estimate = est, 
      `Std. Error` = se,
      `z value` = tval, 
      `Pr(>|z|)` = 2 * pnorm(abs(tval), lower.tail = FALSE))
   }else{
    # Finite sample Student-t
    ans$coefficients <- cbind(
      Estimate = est, 
      `Std. Error` = se,
      `t value` = tval, 
      `Pr(>|t|)` = 2 * pt(abs(tval), rdf, lower.tail = FALSE))
  }
  ans$rdf <- rdf
  ans$sigSq_g <- object$sigSq_g
  ans$sigSq_e <- object$sigSq_e
  ans$delta <- object$delta
  ans$method <- object$method
  ans$formula <- object$formula
  ans$family <- family(object)
  ans$dispersion <- dispersion(object)

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
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#'
#' terms(fit)
#
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
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#'
#' vcov(fit)
#
#' @export
vcov.fastlmm <- function(object, ...) {
  object$vcov
}





#' Extract Model Residuals
#' 
#' Extract Model Residuals
#' 
#' @param object fitted model of class \code{fastlmm}
#' @param type the type of residuals which should be returned.  
#' @param ... other args, not used
#'
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#'
#' residuals(fit)[1:3]
#
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
residuals.fastglmm <- function(object, type = c("deviance" , "pearson", "working", "response"), ...) {

  type <- match.arg(type)

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
    working = 
      (y - mu) / fam$variance(mu), 
    response = y - mu
    )
}




# Compute dispersion of GLM
# adapted from summary.glm()
#' Overdispersion parameter
#'
#' Overdispersion parameter
#'
#' @param object model fit 
#' 
#' @export
dispersion <- function(object){

  df.r <- df.residual(object)
  fam <- family(object)

  if (!is.null(fam$dispersion) && !is.na(fam$dispersion)){
    disp <- fam$dispersion
  }else if(fam$family %in% c("poisson", "binomial")){
    disp <- 1
  }else{

    if( df.r > 0){
      # multiply by scale of weights
      s <- ifelse(is(object, "fastlmm"), object$w.mean, 1)
      w <- object$weights * s
      r <- residuals(object, "working")
      disp <- sum((w*r^2)[w > 0]) / df.r
    } else {
      disp <- NaN
    }
  }

  disp
}







