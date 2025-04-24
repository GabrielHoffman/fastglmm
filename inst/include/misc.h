#ifndef FASTLMM_MISC_H_
#define FASTLMM_MISC_H_

// if -D ARMA, use plain armadillo library
#ifdef ARMA
#include <armadillo>
#else
#include <RcppArmadillo.h>
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

static vec qnorm( const vec & v, const double &mean=0, const double &sd=1 ){
  vec tmp(v);
  for(int i=0; i<tmp.n_elem; i++){
    tmp[i] = R::qnorm(tmp[i], mean, sd, 1, 0);
    // tmp[i] = glm::qnorm(tmp[i],  mean, sd, 1, 0);
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


#endif