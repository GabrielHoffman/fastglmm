#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// only depends on armadillo 

#ifndef DECOMP_H_
#define DECOMP_H_

using namespace arma;

namespace fastlmmLib {


template <typename T1> 
class decomp {

  public:
    decomp( const T1 &Z_, const vec &weights_);

    T1 get_vectors(){return U;}
    vec get_values(){return s;}

  private:
   T1 U;
   vec s;
};

template <typename T1> 
decomp<T1>::decomp( const T1 &Z_, const vec &weights_){

	// R code
	// if (!is.null(weights)) {
	//     Z.mod = sqrt(weights) * Z.mod
	// }
	// cs = colSums(Z.mod^2)
	// idx = order(cs, decreasing = TRUE)
	// cs = cs[idx]
	// Z.mod = Z.mod[, idx]
	// vectors = Z.mod %*% Diagonal(ncol(Z.mod), 1/sqrt(cs))
	// if (!is.null(rank) && rank < ncol(vectors) && rank > 0) {
	//     U <- vectors[, seq_len(rank), drop = FALSE]
	//     cs <- cs[seq_len(rank), drop = FALSE]
	// }
	// list(vectors = vectors, values = as.numeric(cs))

	// https://stackoverflow.com/questions/26204016/in-armadillo-on-c-sumsp-mat-dim-on-sparse-matrices-does-not-work
	// colSums(Z.mod^2)
	// mat sumRows = ones(1, A.n_rows) * A.t() * A;



}


}




#endif