#ifndef LMM_FIT_RESPONSE_H_
#define LMM_FIT_RESPONSE_H_

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
    lmmFitResponses(const T2 &X_, 
                    const T3 &Z_,
                    const double &left_,
                    const double &right_,
                    const double &tol_,
                    const int &nthreads_);

    ModelFitLMMList eval(
                    const T1 &Y,
                    const vector<string> &ids,
                    const mat &Weights);

    private:
    T2 X;  
    T3 Z;
    double left, right, tol;
    int nthreads;
};



// constructor
template <typename T1, typename T2, typename T3> 
lmmFitResponses<T1, T2, T3>::lmmFitResponses(
                            const T2 &X_, 
                            const T3 &Z_,
                            const double &left_,
                            const double &right_,
                            const double &tol_,
                            const int &nthreads_){
    this->X         = X_;
    this->Z         = Z_;
    this->left      = left_;
    this->right     = right_;
    this->tol       = tol_;
    this->nthreads  = nthreads_;
}



template <typename T1, typename T2, typename T3> 
ModelFitLMMList 
  lmmFitResponses<T1, T2, T3>::eval(
                    const T1 &Y,
                    const vector<string> &ids,
                    const mat &Weights){

    // store results
    ModelFitLMMList result(Y.n_cols, ModelFitLMM());

    for( int j = 0; j < Y.n_cols; j++){

        T1 y = Y.col(j);
        vec w = Weights.col(j);

        spectralDecomp<T3> dcmp;
        dcmp.initWithIndicator(Z, w);

        fastlmm fit = fastlmm<T1, T2, T3>(y, X, dcmp.get_vectors(), dcmp.get_values(), w);

        fit.estimate_delta( left, right, tol );

        result.at(j) = fit.get_result();
        result.at(j).ID = ids[j];
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