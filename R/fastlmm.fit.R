# See Design docs:
# https://docs.google.com/document/d/1xO0Gdb5mlgP_54uWmpB_iilabDxAf_OQNXfz7J0EWc0/edit

# Option to pass in pre-computed:
# Xu, Yu, cp_X_low, cp_X_low_Y_low

# Add argument for prior weights?

# Approximate df and rdf
#  just add fixed effect rank to random effect rank
#  assumes they are orthogonal

# eventually, call this _fastlmm
# and call it from within the PQL iteration

# can I estimate theta within PQL?

#' fastglmm
#'
#' Efficient Fit of Generalized Linear Mixed Model with a Single Random Effect
#'
#' @rdname fastglmm-pkg
#' @name fastglmm
#' @useDynLib fastglmm
#' @importFrom Rcpp evalCpp
NULL

#' Convert list to fastlmm class
#'
#' Convert list to fastlmm class
#'
#' @param x list from \code{.fastlmm_()}
#' @param design design matrix for fixed effects
#' @param offset offset
#' @param method method used in model fit
#'
#' @return object of class \code{fastlmm}
#' @export
as.fastlmm <- function(x, design, offset, method) {
  # format results
  x$coefficients <- as.numeric(x$coefficients)
  x$se <- as.numeric(x$se)
  x$design <- design
  x$offset <- offset

  names(x$coefficients) <- colnames(x$design)
  names(x$se) <- colnames(x$design)
  rownames(x$vcov) <- colnames(x$design)
  colnames(x$vcov) <- colnames(x$design)
  x$rank <- ncol(x$design)

  # adapt this to be rdf from H
  x$df.residual <- df.residual(x)
  x$method <- method

  # set names for hessian of variance components
  colnames(x$hessian.vc) <- c("sigSq_g", "sigSq_e")
  rownames(x$hessian.vc) <- c("sigSq_g", "sigSq_e")
  
  class(x) <- "fastlmm"
  x
}


#' Fitter Function for Linear Mixed Model
#'
#' Prepare data for model fitting with a call to Rcpp code
#'
#' @param y response vector
#' @param X design matrix
#' @param Z sparse matrix of indicators for random effect
#' @param offset offset
#' @param rank rank of random effect.  The maximum rank is the number of columns in \code{Z}.  A low rank approximation can be useful if the eigen-values decrease quickly.
#' @param weights an optional vector of prior weights with a value for each sample.  When the response has multiple columns, a vector of weight can be reused for each respose, or a matrix the same dimension as the responses matrix can weight each response separately.
#' @param REML logical scalar - Should the estimates be chosen to optimize the REML criterion vs ML?
#' @param delta  if \code{NULL} estimate delta, if value is given used this fixed values
#' @param delta.range min and max values (in log space), of the search space for delta to fit the random effect
#' @param tol convergence criterion for the 1D search of the delta space
#' @param lambda ridge shrinkage parameter
#' @param nthreads number of threads
#
#' @details Fit a linear mixed model with a single variance component.
#'
#' @return \code{U} and \code{s} values are from the SVD of weighted Z
#'
# other args: sig_a_fixed = FALSE
#' @importFrom methods is
#' @export
fastlmm.fit <- function(y, X, Z, offset = NULL, REML = FALSE, delta = NULL, rank = ncol(Z), weights = NULL, delta.range = c(-10, 10), tol = 1e-6, lambda = 0, nthreads = 6) {

  if (delta.range[1] >= delta.range[2]) {
    stop("delta.range are not valid")
  }

  stopifnot( tol > 0)

  delta.range <- as.numeric(delta.range)

  # add data checks here
  if( ! is.numeric(y) ){
    stop("Response must be numeric")
  }
  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }
 
  if (length(y) != nrow(X)) {
    stop("dimension of Y and X do not match")
  }

  # if delta is NULL, estimate its value
  # by setting to -1 for C++ call
  delta <- ifelse(is.null(delta), -1, delta)

  if (!identical(length(y), length(weights))) {
    stop("Dimension of Y and weights must be the same")
  }

  if (!is.null(offset)) {
    offset <- as.matrix(offset)
    y <- y - offset
  }

  # apply weights
  dcmp <- indicator_decomp(Z, rank = rank, sort=FALSE)

  if (is(Z, "sparseMatrix")) {
    fxn <- .fastlmm_ms
  } else {
    fxn <- .fastlmm_mm
  }

  res <- fxn(
    y = y,
    X = X,
    U = dcmp$vectors,
    s = dcmp$values,
    weights = weights,
    dcmpMethod = "categorical",
    REML = REML,
    delta = delta,
    left = delta.range[1],
    right = delta.range[2],
    tol = tol,
    lambda = lambda,
    nthreads = nthreads
  )  

  res <- as.fastlmm(res, 
          design = X,
          offset = offset, 
          method = ifelse(REML, "REML", "ML"))

  res$Z <- Z

  res
}
