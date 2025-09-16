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

  class(x) <- "fastlmm"
  x
}


#' @export
print.fastlmmList <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(attr(x, "call")), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("Responses:\n")
  coolcat("names(%d): %s\n", names(x))
}



#' Fitter Function for Linear Mixed Model
#'
#' Prepare data for model fitting with a call to Rcpp code
#'
#' @param Y response vector, or matrix with responses as _columns_
#' @param X design matrix
#' @param Z sparse matrix of indicators for random effect
#' @param offset offset
#' @param rank rank of random effect.  The maximum rank is the number of columns in \code{Z}.  A low rank approximation can be useful if the eigen-values decrease quickly.
#' @param weights an optional vector of prior weights with a value for each sample.  When the response has multiple columns, a vector of weight can be reused for each respose, or a matrix the same dimension as the responses matrix can weight each response separately.
#' @param REML logical scalar - Should the estimates be chosen to optimize the REML criterion vs ML?
#' @param delta  if \code{NULL} estimate delta, if value is given used this fixed values
#' @param delta.range min and max values (in log space), of the search space for delta to fit the random effect
#' @param tol convergence criterion for the 1D search of the delta space
#' @param nthreads number of threads
#
#' @details Fit a linear mixed model with a single variance component.
#'
#' @return fill in
#'
# other args: sig_a_fixed = FALSE
#' @importFrom methods is
#' @export
fastlmm.fit <- function(Y, X, Z, offset = NULL, REML = FALSE, delta = NULL, rank = ncol(Z), weights = NULL, delta.range = c(-10, 10), tol = 1e-6, nthreads = 6) {

  if (delta.range[1] >= delta.range[2]) {
    stop("delta.range are not valid")
  }

  stopifnot( tol > 0)

  delta.range <- as.numeric(delta.range)

  # add data checks here
  if (!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }
  if (is.null(weights)) {
    weights <- matrix(1, nrow(Y), ncol(Y))
  }
  if (!is.matrix(weights)) {
    weights <- as.matrix(weights)
  }

  if (nrow(Y) != nrow(X)) {
    stop("dimension of Y and X do not match")
  }

  # if delta is NULL, estimate its value
  # by setting to -1 for C++ call
  delta <- ifelse(is.null(delta), -1, delta)

  if (!identical(dim(Y), dim(weights))) {
    stop("Dimension of Y and weights must be the same")
  }

  if (!is.null(offset)) {
    offset <- as.matrix(offset)
    Y <- Y - offset
  }

  # if 1 response
  if (ncol(Y) == 1) {
    # apply weights
    dcmp <- indicator_decomp(Z, c(weights), rank)

    if (is(Z, "sparseMatrix")) {
      fxn <- .fastlmm_vms
    } else {
      fxn <- .fastlmm_vmm
    }

    res <- fxn(
      y = Y,
      X = X,
      U = dcmp$vectors,
      s = dcmp$values,
      weights = weights,
      REML = REML,
      delta = delta,
      left = delta.range[1],
      right = delta.range[2],
      tol = tol,
      nthreads = nthreads
    )  

    res <- as.fastlmm(res, design = X, offset = offset, method = ifelse(REML, "REML", "ML"))

    # include indicator matrix and its decomposition
    # Does this need to stay in?
    res$Z <- Z
    res$U <- dcmp$vectors
    res$s <- dcmp$values
  } else {
    if (length(weights) == nrow(Y)) {
      weights <- matrix(weights, nrow = nrow(Y), ncol = ncol(Y))
    }

    if( is.null(colnames(Y)) ){
      colnames(Y) = paste0("response_", seq(ncol(Y)))
    }

    if (is(Z, "sparseMatrix")) {
      fxn <- .fastlmm_mms
    } else {
      fxn <- .fastlmm_mmm
    }

    res <- fxn(
      Y = Y,
      ids = colnames(Y),
      X = X,
      Z = Z,
      Weights = weights,
      REML = REML,
      left = delta.range[1],
      right = delta.range[2],
      tol = tol,
      nthreads = nthreads
    )

    # convert each entry to an fastlmm object
    res <- lapply(res, as.fastlmm, design = X, offset = offset, method = ifelse(REML, "REML", "ML"))
    names(res) <- colnames(Y)
    class(res) <- "fastlmmList"
  }

  res
}
