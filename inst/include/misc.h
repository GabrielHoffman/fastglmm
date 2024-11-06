#ifndef FASTLMM_MISC_H_
#define FASTLMM_MISC_H_

// if -D ARMA, use plain armadillo library
#ifdef ARMA
#include <armadillo>
#else
#include <RcppArmadillo.h>
#endif


using namespace arma;

namespace fastlmmLib {


// for each column, scale by w
inline mat scaleRows(const mat &X, const vec &w){
  return X.each_col() % w;  
}


// for each column, scale by w
inline sp_mat scaleRows(const sp_mat &X, const vec &w){

  sp_mat M = sp_mat(X);
  for(size_t i=0; i<X.n_cols; i++){
    M.col(i) %= w;
  }
  return( M );
}


// template <typename T> 
// T scaleRows(const T &X, const vec &w){

//   T M = T(X);
//   for(size_t i=0; i<X.n_cols; i++){
//     M.col(i) %= w;
//   }
//   return( M );
// }

inline mat scaleCols(const mat &X, const vec &w){
  return X.each_row() % w;  
}


// for each column, scale by w
inline sp_mat scaleCols(const sp_mat &X, const vec &w){

  sp_mat M = sp_mat(X);
  for(size_t i=0; i<X.n_cols; i++){
    M.col(i) *= w(i);
  }
  return( M );
}

// template <typename T> 
// T scaleCols(const T &X, const vec &w){

//   T M = T(X);
//   for(size_t i=0; i<X.n_cols; i++){
//     M.col(i) *= w(i);
//   }
//   return( M );
// }

// // scale by w and s
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


// general case, returns false
template <class T>
inline bool isSpMatrix(const T &t) { return false;  } 

 // but for sp_mat returns true
template <>
inline bool isSpMatrix( const sp_mat &t) { return true; } 




}


#endif