#' Spectral decomposition of factor indicator matrix
#'
#' Given a factor, construct the spectral decompostion of the corresponding indicator matrix.  Uses linear time algorithm and sparse matrix algebra.  Resulting vector space is also sparse.
#'
#' @param x object of type \code{factor} or \code{sparseMatrix}
#' @param weights vector of weights with a value for each sample.  If ommited, weights are set to 1.
#' @param rank rank of random effect.  The maximum rank is the number of columns in \code{Z}.  A low rank approximation can be useful if the eigen-values decrease quickly.
#' @param sort sort eigen vectors and values
#'
#' @details This approach is dramatically faster than the naive algorithm that is quadratic time in the number of levels.
#'
#' @return list storing spectral decomposition like from \code{eigen()}
#'
#' @examples
#' set.seed(2)
#'
#' # get 10 samples from 4 individuals
#' indiv <- factor(sample(LETTERS[seq(4)], 10, replace = TRUE))
#'
#' # fast spectral decomposition
#' indicator_decomp(indiv, sort = TRUE)
#'
#' # Create spectral decomposition
#' # using a slower approach ignoring the
#' # special structure of the problem.
#' # This is quadratic time in the number of levels
#' library(Matrix)
#' Z <- t(fac2sparse(indiv))
#' ZtZ <- crossprod(Z)
#' dcmp <- eigen(ZtZ, symmetric = TRUE)
#' dcmp$vectors <- Matrix(zapsmall(dcmp$vectors), sparse = TRUE)
#' A <- solve(dcmp$vectors)
#' D <- Diagonal(length(dcmp$values), 1 / sqrt(dcmp$values))
#' dcmp$vectors <- tcrossprod(Z, A) %*% D
#'
#' dcmp
#' @importFrom Matrix Diagonal t colSums
#' @seealso \code{preprocess_indicator()}
#' @export
indicator_decomp <- function(x, weights = NULL, rank = NULL, sort = FALSE) {
  if (is.factor(x)) {
    Z <- preprocess_indicator(x)
  } else if (is(x, "sparseMatrix")) {
    Z <- x
  }

  if (is.null(weights)) {
    weights <- rep(1, nrow(Z))
  }

  # Compute eigen values and vectors
  evalues <- as.numeric(weights %*% Z)
  names(evalues) <- colnames(Z)
  D <- Diagonal(length(evalues), 1 / sqrt(evalues), names=TRUE)
  vectors <- (Z * sqrt(weights)) %*% D

  if (sort) {
    # sort by cs value
    idx <- order(evalues, decreasing = TRUE)
    evalues <- evalues[idx]
    vectors <- vectors[, idx]
  }

  if (!is.null(rank) && rank < ncol(vectors) && rank > 0) {
    idx <- seq_len(rank)
    U <- vectors[, idx, drop = FALSE]
    evalues <- evalues[idx, drop = FALSE]
  }

  list(vectors = vectors, values = evalues)
}




#' Create sparse indicator matrix
#'
#' Create sparse indicator matrix from factor
#'
#' @param x a \code{factor}
#'
#' @return sparse indicator with levels as columns
#' @importFrom Matrix fac2sparse
#' @export
preprocess_indicator <- function(x) {
  stopifnot(is.factor(x))

  x <- droplevels(x)
  Z.mod <- t(fac2sparse(x))

  Z.mod
}


#' Reconstruct indicator matrix from eigen decomp
#'
#' Reconstruct indicator matrix from eigen decomp
#'
#' @param dcmp eigen decomp from \code{indicator_decomp()}
#' @param weights vector of weights with a value for each sample.  If ommited, weights are set to 1.
#'
#' @examples
#' ID <- factor(sample(LETTERS[1:4], 100, replace=TRUE))
#' w <- seq(length(ID))
#' Z <- preprocess_indicator(ID)
#' dcmp <- indicator_decomp(ID, w)
#' 
#' Z_recon <- reconstruct_indicator( dcmp, w)
#' 
#' range(Z_recon - Z)
#
#' @export
#' @keywords internal
reconstruct_indicator = function( dcmp, weights = NULL ){

  if( is.null(weights) ){
    weights <- rep(1, nrow(dcmp$vectors))
  }

  Z <- with(dcmp, Diagonal(length(weights), 1/sqrt(weights)) %*% vectors %*% Diagonal(length(values), sqrt(values)))

  colnames(Z) <- colnames(dcmp$vectors)

  Z
}





