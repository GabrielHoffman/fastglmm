#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// only depends on armadillo 

#ifndef SPECTRAL_DECOMP_H_
#define SPECTRAL_DECOMP_H_

#include "misc.h"

using namespace arma;

namespace fastlmmLib {


template <typename T> 
class spectralDecomp {

  public:
    spectralDecomp( const T &Z_, const vec &weights_);

    T get_vectors(){return U;}
    vec get_values(){return s;}

  private:
   T U;
   vec s;
};

template <typename T> 
spectralDecomp<T>::spectralDecomp( const T &Z_, const vec &weights_){

	// R code
	// Z.mod = sqrt(weights) * Z.mod
	// cs = colSums(Z.mod^2)
	// vectors = Z.mod %*% Diagonal(ncol(Z.mod), 1/sqrt(cs))

	// list(vectors = vectors, values = as.numeric(cs))

	s = (weights_.t() * Z_).t();
	U = scaleRowsCols(Z_, sqrt(weights_), 1 / sqrt(s));

	// Original
	// T Zmod = scaleRows(Z_, sqrt(weights_));
	// s = Zmod.t() * Zmod * ones(Zmod.n_cols, 1);
	// sp_mat I = speye(Z_.n_cols, Z_.n_cols);
	// U = Zmod * scaleRows(I, 1/sqrt(s));
	
	// s.brief_print();
	// U.brief_print();

	// mat One = mat(1, Z_.n_cols, fill::ones);
	// mat M1 = One * Z_.t();
	// sp_mat I2 = speye(weights_.size(), weights_.size());
	// sp_mat W = scaleRows(I2, weights_);
	// s = (M1 * W * Z_).t();
	// U = scaleRows(Z_, sqrt(weights_)) * scaleRows(I, 1/sqrt(s));


	// mat One = mat(1, Z_.n_cols, fill::ones);
	// mat M1 = One * Z_.t();
	// mat M1 = mat(1, Z_.n_rows, fill::ones)
	// s = ((M1 % weights_.t()) * Z_).t();
	// s = (weights_.t() * Z_).t();
	// sp_mat I = speye(Z_.n_cols, Z_.n_cols);
	// U = scaleRows(Z_, sqrt(weights_)) * scaleRows(I, 1/sqrt(s));
	// U = scaleRowsCols(Z_, sqrt(weights_), 1 / sqrt(s));

	// for( int i=0; i<100; i++){
	// 	// mat One = mat(1, Z_.n_cols, fill::ones);
	// 	// mat M1 = One * Z_.t();
	// 	// s = ((M1 % weights_.t()) * Z_).t();
	// 	// U = scaleRows(Z_, sqrt(weights_)) * scaleRows(I, 1/sqrt(s));
	// 	U = scaleRowsCols(Z_, sqrt(weights_), 1 / sqrt(s));

	// }


}


}




#endif