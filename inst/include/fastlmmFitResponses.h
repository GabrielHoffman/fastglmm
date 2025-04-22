#ifndef FASTLMM_BATCH_RESPONSE_H_
#define FASTLMM_BATCH_RESPONSE_H_

using namespace arma;
using namespace std;

#include "fastlmm.h"

#include "spectralDecomp.h"

namespace fastlmmLib {

// Order of template variables
// T1 Y
// T2 X
// T3 Z
template <typename T1, typename T2, typename T3> 
class lmmFitResponses {

    public:
    lmmFitResponses(const T1 &Y_all_, 
                    const T2 &X_, 
                    const T3 &Z_,
                    const mat &Weights_,
                    const double &left_,
                    const double &right_,
                    const double &tol_,
                    const int &nthreads_);

    ModelFitLMMList eval();

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
lmmFitResponses<T1, T2, T3>::lmmFitResponses(
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



// NOTE: Do not use Rcpp in parallel section
// "C stack usage is too close to the limit"
template <typename T1, typename T2, typename T3> 
ModelFitLMMList 
  lmmFitResponses<T1, T2, T3>::eval(){

    int n_responses = Y_all.n_cols;

    // store results
    ModelFitLMMList result(n_responses, ModelFitLMM());

    for( int i = 0; i < n_responses; i++){

        T1 y = Y_all.col(i);
        vec w = Weights.col(i);

        spectralDecomp<T3> dcmp;
        dcmp.initWithIndicator(Z, w);

        fastlmm fit = fastlmm<T1, T2, T3>(y, X, dcmp.get_vectors(), dcmp.get_values(), w);

        fit.estimate_delta( left, right, tol );

        // #pragma omp critical
        result.at(i) = fit.get_result();
    }
  
    return result;
}


// fastlmm
// Yw = Y_all.col(i) * Weights.col(i),
// Xw = X_orig * Weights.col(i),
// [U, s] = indicator_decomp( Z , Weights.col(i))
// Yu = U_.t() * Yw;
// Xu = U_.t() * Xw;
// cp_X_low = Xw.t() * Xw - Xu.t() * Xu;
// cp_X_low_Y_low = Xw.t() * Yw - Xu.t() * Yu;
// inv_s_delta_Xu = mat( Xu.n_rows, Xu.n_cols);

}




#endif