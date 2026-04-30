#ifndef FASTLMM_MISC_H_
#define FASTLMM_MISC_H_
#pragma once

// if -D USE_R, use RcppArmadillo library
#ifdef USE_R
// [[Rcpp::depends(RcppParallel)]]  
#include <RcppArmadillo.h>
#else
#include <armadillo>
#endif

using namespace arma;


// for each column, scale by w
inline mat scaleEachCol(const mat &X, const vec &w){
  return X.each_col() % w;  
}

// for each column, scale by w
inline sp_mat scaleEachCol(const sp_mat &X, const vec &w){

  sp_mat M = sp_mat(X);
  for(size_t i=0; i<X.n_cols; i++){
    M.col(i) %= w;
  }
  return( M );
}

// for each row, scale by w
inline mat scaleEachRow(const mat &X, const vec &w){
  return w.t() % X.each_row();  
}

// for each row, scale by w
inline sp_mat scaleEachRow(const sp_mat &X, const vec &w){

  sp_mat M = sp_mat(X);
  for(size_t i=0; i<X.n_cols; i++){
    M.col(i) *= w[i];
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


// general case, returns false
template <class T>
inline bool isSpMatrix(const T &t) { return false;  } 

 // but for sp_mat returns true
template <>
inline bool isSpMatrix( const sp_mat &t) { return true; } 


// Given matrix, set entries in rows idx to zero
static void rows_to_zero( mat &X, const uvec &idx ){
  X.rows(idx).zeros();
}

static void rows_to_zero( sp_mat &X, const uvec &idx ){

  for(auto j: idx ){
    for (auto  it = X.begin_row(j); it != X.end_row(j); ++it) {
      *it = 0.0;
    }
  }
  X.clean(0.0);
}


// Return indices (0-based) of rows that contain at least one NaN
static uvec rows_with_nan(const mat& M) {
  std::vector<uword> nan_rows;

  for (uword i = 0; i < M.n_rows; ++i) {
    bool has_nan = false;

    // iterate over elements in row i
    for (mat::const_row_iterator it = M.begin_row(i);
         it != M.end_row(i); ++it) {
      if (std::isnan(*it)) {
        has_nan = true;
        break;  // no need to check rest of row
      }
    }

    if (has_nan){
      nan_rows.push_back(i);
    }
  }

  return conv_to<uvec>::from(nan_rows);
}


// Return indices (0-based) of rows in sp_mat containing NaN
static uvec rows_with_nan(const sp_mat& M) {

  std::vector<uword> nan_rows;

  for (uword i = 0; i < M.n_rows; ++i) {
    bool has_nan = false;

    // iterate over stored elements in row i
    for (sp_mat::const_row_iterator it = M.begin_row(i);
    it != M.end_row(i); ++it){
      if (std::isnan(*it)) {
        has_nan = true;
        break;  // stop once a NaN is found
      }
    }

    if (has_nan)
      nan_rows.push_back(i);
    }

  return conv_to<uvec>::from(nan_rows);
}


// template <typename T>
// const T& min(const T& a, const T& b) {
//     return (a < b) ? a : b;
// }

// template <typename T>
// const T& max(const T& a, const T& b) {
//     return (a > b) ? a : b;
// }

static vec pmax( const vec &v, const double &value){

  vec tmp(v);
  for(int i=0; i<tmp.n_elem; i++){
    tmp[i] = std::max(tmp[i], value);
  }

  return tmp;
}

static vec pmin( const vec &v, const double &value){

  vec tmp(v);
  for(int i=0; i<tmp.n_elem; i++){
    tmp[i] = std::min(tmp[i], value);
  }

  return tmp;
}

#include "qnorm.h"
static vec qnorm( const vec & v, const double &mean=0, const double &sd=1 ){
  vec tmp(v);
  for(int i=0; i<tmp.n_elem; i++){
    // R implementation of qnorm
    // tmp[i] = R::qnorm(tmp[i], mean, sd, 1, 0);

    // C++ without R dependency
    tmp[i] = qnorm_as241(v[i], mean, sd, 1, 0);
  }

  return tmp;
}

static vec y_log_y(const vec & y, const vec & mu){
    // (y) ? (y * log(y/mu)) : 0;

    // initialize to zeros
  vec ret(y.n_elem, fill::zeros);

  for(int i=0; i<y.n_elem; i++){
    if( y[i] != 0.0){
      ret[i] = y[i] * log(y[i]/mu[i]);
    }
  }

  return ret;
}


static double wvar(const vec & x, const vec & w) {

  if (x.n_elem != w.n_elem) {
    throw invalid_argument("Data and weights vectors must have the same number of elements.");
  }

  if (x.is_empty()) {
    return 0.0; 
  }

  double sum_weights = sum(w);

  if (sum_weights == 0) {
    throw runtime_error("Sum of weights cannot be zero.");
  }

  // Calculate weighted mean
  double weighted_mean = sum(x % w) / sum_weights;

  // Calculate weighted variance
  vec diff_sq = square(x - weighted_mean);
  double numerator = sum(w % diff_sq);

  double denominator = sum_weights - (sum(square(w)) / sum_weights);

  if (denominator == 0) {
    return 0.0; 
  }

  return numerator / denominator;
}


// Adapted from tbb::blocked_range
// Designed to be used when tbb is not available
template<typename T>
class blocked_range {
  public:
    // constructors
    blocked_range( T begin, T end ) : 
      begin_value(begin), end_value(end) {}

    // get range limits
    T begin(){ return begin_value; }
    T end(){ return end_value; }
  private:
    T begin_value, end_value;
};


// Set OpenMP threads to 1 and disable nested parallelism
static void disable_parallel_blas(){
  #ifdef _OPENMP
  #include <omp.h> 
  // set threads
  omp_set_num_threads(1);
  // disable nested parallelism
  omp_set_max_active_levels(0);
  #endif
}

// For NAN entries in Y, set matching entry in W to zero
static void match_NAN_zeros( const mat &Y, mat &W){
  for (arma::uword r = 0; r < Y.n_rows; ++r) {
    for (arma::uword c = 0; c < Y.n_cols; ++c) {
      if (std::isnan(Y(r, c))) {
        W(r, c) = 0.0;
      }
    }
  }
}

static int count_nan( const vec &v){
  int nan_count = 0;
  for (arma::uword i = 0; i < v.n_elem; ++i) {
    if (std::isnan(v(i))) {
      nan_count++;
    }
  }
  return nan_count;
}


static uvec countResponseFilter( const mat &Y){

  // # filter genes by expression
  // keep1 <- rowSums2(countMatrix > 0) > max(2, 0.1*ncol(countMatrix))
  // keep2 <- rowSums2(countMatrix > 1) > max(2, 0.01*ncol(countMatrix))
  // keep <- keep1 & keep2

  int nr = Y.n_rows;

  uvec keep1 = find(sum(Y > 0, 0) > max(2, (int) 0.1*nr));
  uvec keep2 = find(sum(Y > 1, 0) > max(2, (int) 0.01*nr));

  return intersect(keep1, keep2);
}

/** return subset of v defined by idx
 */ 
template <typename T>
static vector<T> subset( const vector<T> &v, const uvec &idx){

  vector<T> res;
  res.reserve(idx.n_elem);

  for (uword i = 0; i < idx.n_elem; ++i) {
    res.push_back( v[idx[i]] );
  }

  return res;
}

/** Robust mean using only samples less than z-cuttoff
*/
static double robust_mean( const vec &x, const double &z_cutoff){

  // keep only finite values
  vec x1 = x.elem(find_finite(x));

  // if no elements retained
  if( x1.n_elem == 0){
    return datum::nan;
  }

  // compute z-score
  // add 1e-15 to avoid issue with sd is zero
  vec z = (x1 - mean(x1)) / (stddev(x1) + 1e-15);

  // find indeces where abs z-score is less than cutoff
  uvec idx = find(abs(z) < z_cutoff);

  // mean of retained values
  return mean(x1.elem(idx));

}

/** Check that all values in x are integer values
 */
static bool all_integer_valued(const vec &x, double tol = 1e-12) {
  return all(abs(x - round(x)) < tol);
}



#endif