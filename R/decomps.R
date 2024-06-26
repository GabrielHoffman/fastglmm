


#' Spetral decomposition of factor indicator matrix
#' 
#' Given a factor, construct the spectral decompostion of the corresponding indicator matrix.  Uses linear time algorithm and sparse matrix algebra.  Resulting vector space is also sparse.
#' 
#' @param Factor object of type \code{factor}
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
#' @importFrom Matrix fac2sparse Diagonal t
#' @export 
indicator_decomp = function( Factor ){
	Factor = droplevels( Factor )
	tab = table(Factor)
	tab = tab[order(tab, decreasing=TRUE)]
	lvls = names(tab[order(tab, decreasing=TRUE)])
	ID = factor(Factor, lvls)
	Z.mod = t(fac2sparse(ID, to='l'))
	vectors = Z.mod %*% Diagonal(length(tab), 1/sqrt(as.numeric(tab)))
	list(vectors = vectors, values = as.numeric(tab))
}

