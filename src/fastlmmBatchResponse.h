#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "fastlmm.h"
#include "spectralDecomp.h"

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
// T3 Z
template <typename T1, typename T2, typename T3> 
class fastlmmBatchResponse {

    public:
    fastlmmBatchResponse(   const T1 &Y_all_, 
                            const T2 &X_, 
                            const T3 &Z_,
                            const mat &Weights_,
                            const double &left_,
                            const double &right_,
                            const double &tol_,
                            const int &nthreads_);

    vector<fastlmm_result> eval();

    private:
    T1 Y_all; 
    T2 X;  
    T3 Z;
    mat Weights;
    double left, right, tol;
    int nthreads;
};



// constructor
template <typename T1, typename T2, typename T3> 
fastlmmBatchResponse<T1, T2, T3>::fastlmmBatchResponse(
                            const T1 &Y_all_, 
                            const T2 &X_, 
                            const T3 &Z_,
                            const mat &Weights_,
                            const double &left_,
                            const double &right_,
                            const double &tol_,
                            const int &nthreads_){
    this->Y_all     = Y_all_;
    this->X         = X_;
    this->Z         = Z_;
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
    // fastlmm fit = fastlmm<T1, T2, T3>();

    // iterate through responses 
    #pragma omp for 
    for( int i = 0; i < n_responses; i++){

        spectralDecomp dcmp = spectralDecomp<T3>(Z, Weights.col(i));
        T1 Yw = Y_all.col(i) * Weights.col(i);
        T2 Xw = X % Weights.col(i);        
        T1 Yu = dcmp.get_vectors().t() * Yw;
        mat Xu = dcmp.get_vectors().t() * conv_to<mat>::from(Xw); 
        mat cp_X_low = Xw.t() * Xw - Xu.t() * Xu;
        mat cp_X_low_Y_low = Xw.t() * Yw - Xu.t() * Yu;  


        fastlmm fit = fastlmm(Yw, Xu, dcmp.get_vectors(), dcmp.get_values(), Weights.col(i), Yu, Xu, cp_X_low, cp_X_low_Y_low);

        fit.estimate_delta( left, right, tol );

        #pragma omp critical
        result.at(i) = fit.get_result();
    }
  }

  return result;
}







}




#endif