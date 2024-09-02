#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifndef MISC_H_
#define MISC_H_

using namespace arma;

namespace fastlmmLib {

// for each column, scale by w
mat scaleRows(const mat &X, const vec &w){
  return X.each_col() % w;  
}


// for each column, scale by w
sp_mat scaleRows(const sp_mat &X, const vec &w){

  sp_mat M = sp_mat(X);
  for(size_t i=0; i<X.n_cols; i++){
    M.col(i) %= w;
  }
  return( M );
}



mat scaleCols(const mat &X, const vec &w){
  return X.each_row() % w;  
}


// for each column, scale by w
sp_mat scaleCols(const sp_mat &X, const vec &w){

  sp_mat M = sp_mat(X);
  for(size_t i=0; i<X.n_cols; i++){
    M.col(i) *= w(i);
  }
  return( M );
}

// scale by w and s
template <typename T> 
T scaleRowsCols(const T &X, const vec &w1, const vec &w2){

  T M = T(X);
  for(size_t i=0; i<X.n_cols; i++){
    // M.col(i) %= w1 * w2(i);
    // manually setting order of operations is _much_ faster
    // since this reduces the number of multiplications
    M.col(i) *= w2(i);
    M.col(i) %= w1;
  }
  return( M );
}




}


#endif