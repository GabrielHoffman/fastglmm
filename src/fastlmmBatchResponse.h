#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "fastlmm.h"

#ifdef _OPENMP
    // [[Rcpp::plugins(openmp)]]
    #include <omp.h>
#else
    #define omp_get_num_threads() 0
    #define omp_get_thread_num() 0
#endif


#ifndef FASTLMM_BATCH_RESPONSE_H_
#define FASTLMM_BATCH_RESPONSE_H_

using namespace arma;

namespace fastlmmLib {

// Order of template variables
// T1 Y
// T2 X
// T3 U
template <typename T1, typename T2, typename T3> 
class fastlmmBatchResponse {

    public:
    fastlmmBatchResponse(   const T1 &Y_all_, 
                            const T2 &X_, 
                            const T3 &U_, 
                            const vec &s_,
                            const mat &Weights_,
                            const double &left_,
                            const double &right_,
                            const double &tol_,
                            const int &nthreads_);

    vector<fastlmm_result> eval();

    private:
    T1 Y_all; 
    T2 X;  
    T3 U; 
    vec s;
    mat Weights;
    double left, right, tol;
    int nthreads;
};



// constructor
template <typename T1, typename T2, typename T3> 
fastlmmBatchResponse<T1, T2, T3>::fastlmmBatchResponse(
                            const T1 &Y_all_, 
                            const T2 &X_, 
                            const T3 &U_, 
                            const vec &s_,
                            const mat &Weights_,
                            const double &left_,
                            const double &right_,
                            const double &tol_,
                            const int &nthreads_){
    this->Y_all     = Y_all_;
    this->X         = X_;
    this->U         = U_;
    this->s         = s_;
    this->Weights   = Weights_;
    this->left      = left_;
    this->right     = right_;
    this->tol       = tol_;
    this->nthreads  = nthreads_;
}


template <typename T1, typename T2, typename T3> 
vector<fastlmm_result> 
  fastlmmBatchResponse<T1, T2, T3>::eval(){

  Rcpp::Rcout << "Fit batch response" << std::endl;

  // need to apply weights matrix Y_all_, decomp, and X
  mat Yu_all = U.t() * Y_all;
  int n_responses = Y_all.n_cols;

  // store results
  vector<fastlmm_result> result(n_responses, fastlmm_result());

  // NOTE: Do not use Rcpp in parallel section
  // "C stack usage is too close to the limit"

  #ifdef _OPENMP 
  // set threads
  omp_set_num_threads(nthreads);
  // disable nested parallelism
  omp_set_max_active_levels(1);
  #endif

  // fastlmm
  // Yw = Y_all.col(i) * Weights.col(i),
  // Xw = X_orig * Weights.col(i),
  // [U, s] = indicator_decomp( Z , Weights.col(i))
  // Yu = U_.t() * Yw;
  // Xu = U_.t() * Xw;
  // cp_X_low = Xw.t() * Xw - Xu.t() * Xu;
  // cp_X_low_Y_low = Xw.t() * Yw - Xu.t() * Yu;
  // inv_s_delta_Xu = mat( Xu.n_rows, Xu.n_cols);

  #pragma omp parallel
  {
    // initialize
    fastlmm fit = fastlmm<T1, T2, T3>(X, U, s);

    // iterate through responses 
    #pragma omp for 
    for( int i = 0; i < n_responses; i++){

        decomp dcmp = decomp(Z, Weights.col(i));



        // fit.update_response(Y_all.col(i), 
        //                   Weights.col(i), 
        //                   Yu_all.col(i));


        fit.estimate_delta( left, right, tol );

        #pragma omp critical
        result.at(i) = fit.get_result();
    }
  }

  return result;
}







}




#endif