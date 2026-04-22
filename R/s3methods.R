

# Class definitions
###################

#' fastlmm
#'
#' Stores results of fastlmm model fit
#'
#' @name fastlmm-class
#' @rdname fastlmm-class
#' @exportClass fastlmm
#' @return none
#' @keywords internal
setClass("fastlmm", contains="list")


#' fastglmm
#'
#' Stores results of fastglmm model fit
#'
#' @name fastglmm-class
#' @rdname fastglmm-class
#' @exportClass fastglmm
#' @return none
#' @keywords internal
setClass("fastglmm", contains="fastlmm")



#' ANOVA Tables
#' 
#' ANOVA Tables
#' 
#' @param object fitted model of class \code{fastlmm}
#' @param ddf \code{"satterthwaite"}: use Satterthwaite approximation to denominator degrees of freedom for the F distribution, or \code{"asymptotic"} to use chisq distribution as null for the test statistic
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
anova.fastlmm <- function(object, ddf = c("satterthwaite", "asymptotic"), ...){

  ddf <- match.arg(ddf)

  # Assign each coef to a contrast
  asgn <- attr(object$design, "assign")

  # Names of effects
  nmeffects <- attr(terms(object), "term.labels")[unique(asgn)]

  if( attr(terms(object),"intercept") == 1){
    nmeffects <- c("(Intercept)", nmeffects)
  }

  # create list of contrast matricies, one for each variable
  L.list <- lapply(seq(0, max(asgn)), function(i){

    hm <- colnames(object$design)[which(asgn == i)]
    createContrastMatrix( object, hm)
  })

  if( ddf == "satterthwaite"){
    # for each contrast
    df2 <- sapply(L.list, function(L){
      # compute ddf for each coefficient
      nu <- ddf(object, L)

      # if joint contrast, summarize nu values
      get_Fstat_ddf(nu)
        })
  }

  df1 <- Chisq <- NULL

  res <- lapply(seq(0, max(asgn)) + 1, function(i){

    # get contrast matrix 
    L <- L.list[[i]]

    # value being tested
    LBeta <- L %*% coef(object) 

    # its variance
    VLbeta <- L %*% tcrossprod(vcov(object), L)

    # F statistic
    F.stat <- crossprod(LBeta, solve(VLbeta, LBeta)) / ncol(VLbeta)

    if( ddf == "satterthwaite"){
      res <- data.frame(
        id = nmeffects[i], 
        df1 = nrow(L), 
        df2 = df2[i],
        F = F.stat) %>%
      mutate('Pr(>F)' = pf(F, df1, df2, lower.tail=FALSE))
    }else{
      res <- data.frame(
        id = nmeffects[i], 
        df = nrow(L)) %>%
      mutate(
        Chisq = F.stat*df, 
        'Pr(>Chisq)' = pchisq(Chisq, df, lower.tail=FALSE))
    }
    res
  }) %>%
  bind_rows %>%
  column_to_rownames("id") 

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
#' Residual degrees-of-freedom is the trace of the residual hat matrix
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
#' Return diagonals of the hat matrix
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
#' Effective degrees-of-freedom of model fit is the trace of the hat matrix mapping from observed to predicted response
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
#' @param newdata optionally, a data.frame in which to look for variables with which to predict.  If omitted, the fitted  predictors are used.
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
#' @importFrom stats fitted model.frame model.offset
#' @importFrom Matrix crossprod
#' @importFrom reformulas nobars findbars
#' @export
fitted.fastlmm <- function(object, ..., newdata = NULL) {
    
  blp = ranef(object)

  if( is.null(newdata) ){
    # fitted value for existing data
    a <- with(object, U %*% (sqrt(s) * crossprod(V, blp[[1]])))
    v <- as.matrix(a) / sqrt(object$weights) + object$design %*% coef(object)

    if (!is.null(object$offset)) {
      v <- v + object$offset
    }
  }else{
    # fitted value for _new_ data
    design.star <- model.matrix(nobars(object$formula), newdata)

    # random effects matrix
    Z.star <- t(fac2sparse(newdata[[names(blp)]], drop.unused.levels=FALSE))

    ids <- rownames(blp[[1]])

    v <- design.star %*% coef(object) + Z.star[,ids] %*% blp[[1]]

    # offset
    mf <- model.frame(nobars(object$formula), newdata)
    os =  model.offset(mf)
    if (!is.null(os)) {
      v <- v + os
    }
  }
  as.matrix(v)
}


# get eta from fastglmm
get_eta = function(object, newdata = NULL){
  class(object) <- "fastlmm"
  fitted(object, newdata = newdata)
}

#' @rdname fitted
#' @importFrom stats fitted
#' @export
fitted.fastglmm = function(object,..., newdata = NULL){

  # convert eta to mu
  object$family$linkinv( get_eta(object, newdata = newdata) )
}


#' Test Linear Hypothesis
#' 
#' Test Linear Hypothesis
#' 
#' @param model fitted model of class \code{fastlmm}
#' @param hypothesis.matrix matrix (or vector) giving linear combinations of coefficients by rows, or a character vector giving the hypothesis in symbolic form 
#' @param rhs right-hand-side vector for hypothesis, with as many entries as rows in the hypothesis matrix; can be omitted, in which case it defaults to a vector of zeroes. For a multivariate linear model, ‘rhs’ is a matrix, defaulting to 0
#' @param ... other args passed to \code{car::linearHypothesis.default()}
#' @param ddf \code{"satterthwaite"}: use Satterthwaite approximation to denominator degrees of freedom for the Student-t or F distribution, or \code{"asymptotic"} to use normal distribution or chisq as null for the test statistic
#'
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#' 
#' linearHypothesis(fit, "trtdrug")
#
#' @importFrom car linearHypothesis linearHypothesis.default makeHypothesis
#' 
#' @rdname linearHypothesis
#' @seealso \code{car::linearHypothesis()}
#' @export
linearHypothesis <- function (model, ...) 
{
    UseMethod("linearHypothesis")
}

#' @rdname linearHypothesis
#' @export
linearHypothesis.fastlmm <- function(model, hypothesis.matrix, rhs = NULL, ..., ddf = c("satterthwaite", "asymptotic")){

  ddf <- match.arg(ddf)

  L <- createContrastMatrix( model, hypothesis.matrix)

  if( ddf == "satterthwaite" ){
    error.df <- ddf(model, L)
    error.df <- get_Fstat_ddf(error.df)
    test <- "F"
  }else{
    # Asymptotic null distribution
    error.df <- Inf
    test <- "Chisq"
  }

  linearHypothesis.default(model, L, rhs, test=test,..., error.df=error.df)
}





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


#' Extracting the Model Frame from a Fit
#' 
#' Extracting the Model Frame from a Fit
#'
#' @param formula regression fit
#' @param ... other args
#' 
#' @examples
#' library(MASS)
#' 
#' # GLMM via PQL
#' fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
#'    family = binomial(), data = bacteria)
#'
#' data <- model.frame(fit)
#' head(data)
#' @importFrom reformulas subbars
#' @importFrom stats model.frame.default
#' @export
model.frame.fastlmm <- function(formula, ...){

  object <- formula

  X <- model.frame.default(
    subbars(object$formula), 
    object$data,
    ...)

  if( !is.null(object$weights) ){
    X[,'(weights)'] <- object$weights 
  }

  X
}




#' Family function for Negative Binomial GLMs
#' 
#' Specifies the information required to fit a Negative Binomial generalized linear model, with known \code{theta} parameter.  
#'
#' @param theta The known value of the additional parameter, \code{theta}. If set to \code{NA}, \code{fastglmm()} will estimate the parameter from the data
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


#' Diagnostic Plots for Model Fits
#' 
#' Diagnostic plots for \code{fastlmm} Fits
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
#' @param newdata optionally, a data.frame in which to look for variables with which to predict.  If omitted, the fitted  predictors are used.
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

  if (!is.null(newdata) & type == "terms") {
    stop("newdata is not currently supported for type terms")
  }

  switch(type, 
    response = fitted(object, newdata = newdata),
    terms = predict_terms(object))
}


#' @rdname predict
#' @export
predict.fastglmm <- function(object, newdata = NULL,
                 type = c("link", "response", "terms"), ...) {

  type <- match.arg(type)

  if (!is.null(newdata) & type == "terms") {
    stop("newdata is not currently supported for type terms")
  }
  
  switch(type,
      link = get_eta(object, newdata = newdata),
      response = fitted(object, newdata = newdata),
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




#' Extract the modes of the random effect
#'
#' Extract the conditional modes of the random effect from the model fit
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
#' @importFrom Matrix crossprod
#' @method ranef fastlmm
#' @importFrom nlme ranef
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

  lst <- list()
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
#' @method fixef fastlmm
#'
#' @importFrom nlme fixef
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
    cat("Generalized linear mixed model fit by", x$method, "['fastglmm']\n")
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
  printCoefmat(coefs, digits = digits, signif.stars = signif.stars, na.print = "NA", tst.ind = which(colnames(coefs) == "df"), ...)

  # cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ", format(x$dispersion, digits=3), ")\n", sep='')

  cat("\nResidual df:", round(x$rdf, digits=1), "\n")

  cat("\nVariance components:")
  cat("\n  sigSq_g:", format(x$sigSq_g, digits = digits))
  cat("\n  sigSq_e:", format(x$sigSq_e, digits = digits))
  cat("\n  delta:  ", format(x$delta, digits=digits))

  cat("\n")
}




#' Object Summaries and Hypothesis Testing
#'
#' Object summaries and hypothesis testing of fixed effects
#' 
#' @param object fitted model of class \code{fastlmm}
#' @param ddf \code{"satterthwaite"}: use Satterthwaite approximation to denominator degrees of freedom for the Student-t distribution, or \code{"asymptotic"} to use normal distribution as null for the test statistic
#' @param ... other args, not used
#'
#' @importFrom stats coef pt
#' @export
summary.fastlmm <- function(object, ddf = c("satterthwaite", "asymptotic"), ...) {

  ddf <- match.arg(ddf)
  z <- object

  est <- coef(object)
  se <- object$se 
  rdf <- df.residual(object)
  tval <- est / se
  ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
  ans$aliased <- is.na(coef(object))
  ans$residuals <- z$residuals

  if( ddf == "satterthwaite"){
    # Finite sample Student-t, with approximate ddf
    ddf.values <- ddf(object)
    ans$coefficients <- cbind(
      Estimate = est, 
      `Std. Error` = se,
      df = round(ddf.values, 1),
      `t value` = tval, 
      `Pr(>|t|)` = 2 * pt(abs(tval), ddf.values, lower.tail = FALSE))
  }else{
    # Normal approximation
    ans$coefficients <- cbind(
      Estimate = est, 
      `Std. Error` = se,
      `z value` = tval, 
      `Pr(>|z|)` = 2 * pnorm(abs(tval), lower.tail = FALSE))
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
#' @importFrom reformulas nobars
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
#' @rdname vcov
#' @export
vcov.fastlmm <- function(object, ...) {

  object$vcov
}

# #' @rdname vcov
# #' @export
# vcov.fastglmm <- function(object, ...) {

#   # quasi-likelihood dispersion 
#   if( is.na(family(object)$dispersion) ){
#     phi <- object$dispersion
#   }else{
#     phi <- 1
#   }

#   object$vcov * phi
# }



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
 
  os <- 0
  if (!is.null(object$offset)) {
    os <- object$offset
  }

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
      object$y / sqrt(wts) - fam$linkfun(fitted(object)) + os,
    response = y - mu 
    )
}




#' Overdispersion parameter phi for quasi-likelihood
#'
#' Overdispersion parameter phi for quasi-likelihood
#'
#' @param object model fit 
#' 
#' @export
setGeneric("dispersion", 
  function(object) {
  standardGeneric("dispersion")
})


#' @rdname dispersion
#' @export
setMethod("dispersion", signature("fastlmm"), 
  function(object) {

  object$dispersion
})




# #' @rdname dispersion
# #' @export
# setMethod("dispersion", signature("fastlmm"), 
#   function(object) {

#   df.r <- df.residual(object)
#   fam <- family(object)

#   if (!is.null(fam$dispersion) && !is.na(fam$dispersion)){
#     disp <- fam$dispersion
#   }else if(fam$family %in% c("poisson", "binomial")){
#     disp <- 1
#   }else{

#     if( df.r > 0){
#       if(is(object, "fastglmm")){
#         w <- object$prior.weights
#       }else{        
#         # multiply by scale of weights
#         s <- ifelse(is(object, "fastlmm"), object$w.mean, 1)
#         w <- object$weights * s
#       }
#       r <- residuals(object, "pearson")
#       disp <- sum((w*r^2)[w > 0]) / df.r
#     } else {
#       disp <- NaN
#     }
#   }

#   disp
# })







