#ifndef LMM_FIT_FEATURES_H_
#define LMM_FIT_FEATURES_H_

#include "fastlmm.h"
#include "spectralDecomp.h"

using namespace arma;
using namespace std;

namespace fastlmmLib {

// Order of template variables
// T1 Y
// T2 X
// T3 Z
template <typename T1, typename T2, typename T3> 
class lmmFitFeatures {

    public:
    lmmFitFeatures( const T1 &Y_, 
                    const T2 &X_, 
                    const T3 &U_,
                    const vec &s_,
                    const vec &weights_,
                    const double &delta_,
                    const double &left_,
                    const double &right_,
                    const double &tol_,
                    const int &nthreads_);

    ModelFitLMMList eval(const T2 &X_add_,
                        const vector<string> &ids);

    private:
    T1 Y; 
    T2 X_shared;  
    T3 U;
    vec s, weights;
    double delta, left, right, tol;
    int nthreads;
    fastlmm<T1, T2, T3> fit;
    spectralDecomp<T3> dcmp;
};



// constructor
template <typename T1, typename T2, typename T3> 
lmmFitFeatures<T1, T2, T3>::lmmFitFeatures(
                            const T1 &Y_, 
                            const T2 &X_, 
                            const T3 &U_,
                            const vec &s_,
                            const vec &weights_,     
                            const double &delta_,
                            const double &left_,
                            const double &right_,
                            const double &tol_,
                            const int &nthreads_){

    // initialize internal variables
    this->Y         = Y_;
    this->X_shared  = X_;
    // this->U         = U_;
    // this->s         = s_;
    this->weights   = weights_;
    this->delta = delta_;
    this->left = left_;
    this->right = right_;
    this->tol = tol_;
    this->nthreads = nthreads_;

    // curently no reweighting
    // , weights_
    dcmp.initWithEigenDecomp(U_, s_);
}



// NOTE: Do not use Rcpp in parallel section
// "C stack usage is too close to the limit"
template <typename T1, typename T2, typename T3> 
ModelFitLMMList 
  lmmFitFeatures<T1, T2, T3>::eval( const T2 &X_add_,
                                    const vector<string> &ids){

    int n_tests = X_add_.n_cols;

    // store results
    ModelFitLMMList result(n_tests, ModelFitLMM());

    for( int j = 0; j < n_tests; j++){

        // currently, only 1 cbind'd
        mat X_combined = join_horiz(X_shared, X_add_.col(j));

        // fits full model each time,
        // for speed, need to save Y, X, scaled by U and s
        fastlmm fit = fastlmm(Y, X_combined, dcmp.get_vectors(), dcmp.get_values(), weights);

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