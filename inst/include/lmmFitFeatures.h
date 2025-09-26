#ifndef LMM_FIT_FEATURES_H_
#define LMM_FIT_FEATURES_H_

#include "spectralDecomp.h"

using namespace arma;
using namespace std;
using namespace fastglmmLib;

namespace fastglmmLib {

// Order of template variables
// T1 Y
// T2 X
// T3 Z
template <typename T1, typename T2, typename T3> 
class lmmFitFeatures {

  public:
  lmmFitFeatures( const T1 &Y, 
                  const T2 &X, 
                  const T3 &U,
                  const vec &s,
                  const vec &weights,
                  const double &delta,
                  const double &left = -10,
                  const double &right = 10,
                  const double &tol = 1e-6,
                  const int &nthreads = 1,
                  const ModelDetail md = LOW,
                  const bool REML = false);

  ModelFitLMMList eval(const T2 &X_add_,
                        const vector<string> &ids);

  private:
  T1 Y; 
  T2 X_shared;  
  T3 U;
  vec s, weights;
  double delta, left, right, tol;
  int nthreads;
  ModelDetail md;
  bool REML;
  spectralDecomp<T3> dcmp;
};



// constructor
template <typename T1, typename T2, typename T3> 
lmmFitFeatures<T1, T2, T3>::lmmFitFeatures(
                            const T1 &Y, 
                            const T2 &X, 
                            const T3 &U,
                            const vec &s,
                            const vec &weights,     
                            const double &delta,
                            const double &left,
                            const double &right,
                            const double &tol,
                            const int &nthreads,
                            const ModelDetail md,
                            const bool REML) :
  Y(Y), 
  X_shared(X), 
  weights(weights),
  delta(delta), 
  left(left), 
  right(right), 
  tol(tol), 
  nthreads(nthreads), 
  md(md), 
  REML(REML)
  {

  dcmp.initWithEigenDecomp(U, s);
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

  // Parallel part using Thread Building Blocks
  tbb::task_arena limited_arena(nthreads);
  limited_arena.execute([&] {
  tbb::parallel_for(
    tbb::blocked_range<int>(0, n_tests, 100), 
    [&](const tbb::blocked_range<int>& r){ 

    disable_parallel_blas();

    // iterate through responses 
    for (int j = r.begin(); j != r.end(); ++j) { 

      // currently, only 1 cbind'd
      mat X_combined = join_horiz(X_shared, X_add_.col(j));

      // fits full model each time,
      // for speed, need to save Y, X, scaled by U and s
      fastlmm fit = fastlmm(Y, X_combined, dcmp.get_vectors(), dcmp.get_values(), weights, md, REML);

      if( delta > 0 ){
          fit.eval_delta( delta ); 
      }else{
          fit.estimate_delta( left, right, tol );
      }

      result.at(j) = fit.get_result();
      result.at(j).ID = ids[j];
    }
  }); });

  return result;
}


// fastlmm
// Yw = Y_all.col(i) * Weights.col(i),
// Xw = X_orig * Weights.col(i),
// [U, s] = indicator_decomp( Z , Weights.col(i))
// Yu = U_.t() * Yw;
// Xu = U_.t() * Xw;
// Gamma_XX = Xw.t() * Xw - Xu.t() * Xu;
// Gamma_XY = Xw.t() * Yw - Xu.t() * Yu;
// inv_s_delta_Xu = mat( Xu.n_rows, Xu.n_cols);

}




#endif