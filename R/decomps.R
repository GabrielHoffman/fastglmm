


#' Spectral decomposition of factor indicator matrix
#' 
#' Given a factor, construct the spectral decompostion of the corresponding indicator matrix.  Uses linear time algorithm and sparse matrix algebra.  Resulting vector space is also sparse.
#' 
#' @param x object of type \code{factor} or \code{sparseMatrix}
#' @param weights vector of weights with a value for each sample.  If ommited, weights are set to 1.
#' @param rank rank of random effect.  The maximum rank is the number of columns in \code{Z}.  A low rank approximation can be useful if the eigen-values decrease quickly.
#' 
#' @details This approach is dramatically faster than the naive algorithm that is quadratic time in the number of levels.
#' 
#' @return list storing spectral decomposition like from \code{eigen()}
#' 
#' @examples
#' set.seed(2)
#' 
#' # get 10 samples from 4 individuals
#' indiv = factor(sample(LETTERS[seq(4)], 10, replace=TRUE))
#' 
#' # fast spectral decomposition
#' indicator_decomp( indiv )
#' 
#' # Create spectral decomposition
#' # using a slower approach ignoring the 
#' # special structure of the problem.
#' # This is quadratic time in the number of levels
#' library(Matrix)
#' Z = t(fac2sparse(indiv))
#' ZtZ = crossprod(Z)
#' dcmp = eigen(ZtZ, symmetric=TRUE)
#' dcmp$vectors = Matrix(zapsmall(dcmp$vectors), sparse = TRUE)
#' A = solve(dcmp$vectors)
#' D = Diagonal(length(dcmp$values), 1/sqrt(dcmp$values))
#' dcmp$vectors = tcrossprod(Z, A) %*% D
#' 
#' dcmp
#' @importFrom Matrix Diagonal t colSums
#' @export 
indicator_decomp = function( x, weights = NULL, rank = NULL){

	if( is.factor(x) ){
		Z.mod = preprocess_indicator( x )
	}else if( is(x, "sparseMatrix") ){
		Z.mod = x
	}

	# weight the rows of the indicator matrix
	if( ! is.null(weights) ){
		Z.mod = sqrt(weights) * Z.mod
	}

	# compute col sum of squares
	cs = colSums(Z.mod^2)

	# sort by cs value
	idx = order(cs, decreasing=TRUE)
	cs = cs[idx]
	Z.mod = Z.mod[,idx]

	vectors = Z.mod %*% Diagonal(ncol(Z.mod), 1/sqrt(cs))

	if( !is.null(rank) && rank < ncol(vectors) && rank > 0){
		U <- vectors[,seq_len(rank), drop=FALSE]
		cs <- cs[seq_len(rank), drop=FALSE]
	}

	list(vectors = vectors, values = as.numeric(cs))
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
 preprocess_indicator = function( x ){
 	stopifnot( is.factor(x) )

<<<<<<< Updated upstream
#' Create sparse indicator matrix
#' 
#' Create sparse indicator matrix from factor
#'
#' @param x a \code{factor}
#' 
#' @return sparse indicator with levels as columns
#' @importFrom Matrix fac2sparse
#' @export 
preprocess_indicator = function( x ){
	stopifnot( is.factor(x) )

	x = droplevels( x )
	Z.mod = t(fac2sparse(x))

	Z.mod
}
=======
 	x = droplevels( x )
 	Z.mod = t(fac2sparse(x))
>>>>>>> Stashed changes

 	Z.mod
 }









