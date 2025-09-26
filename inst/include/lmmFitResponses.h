#ifndef LMM_FIT_RESPONSE_H_
#define LMM_FIT_RESPONSE_H_

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

#include "fastlmm_fit.h"
#include "ModelFit.h"
#include "spectralDecomp.h"

using namespace arma;
using namespace std;

namespace fastglmmLib {

// Order of template variables
// T1 Y
// T2 X
// T3 Z
template <typename T1, typename T2, typename T3> 
class lmmFitResponses {

  public:
  lmmFitResponses(const T2 &X, 
                  const T3 &Z,
                  const double &left = -10,
                  const double &right = 10,
                  const double &tol = 1e-6,
                  const int &nthreads = 1,
                  const ModelDetail md = LOW,
                  const bool REML = false);

  ModelFitLMMList eval(
                  const T1 &Y,
                  const vector<string> &ids,
                  const mat &Weights);

  private:
  T2 X;  
  T3 Z;
  double left, right, tol;
  int nthreads;
  ModelDetail md;
  bool REML;
};



// constructor
template <typename T1, typename T2, typename T3> 
lmmFitResponses<T1, T2, T3>::lmmFitResponses(
                            const T2 &X, 
                            const T3 &Z,
                            const double &left,
                            const double &right,
                            const double &tol,
                            const int &nthreads,
                            const ModelDetail md,
                            const bool REML):
  X(X), 
  Z(Z),
  left(left),
  right(right),
  tol(tol),
  nthreads(nthreads),
  md(md),
  REML(REML) {
}



template <typename T1, typename T2, typename T3> 
ModelFitLMMList 
  lmmFitResponses<T1, T2, T3>::eval(
                    const T1 &Y,
                    const vector<string> &ids,
                    const mat &Weights){

  // store results
  ModelFitLMMList result(Y.n_cols, ModelFitLMM());

  // Parallel part using Thread Building Blocks
  tbb::task_arena limited_arena(nthreads);
  limited_arena.execute([&] {
  tbb::parallel_for(
    tbb::blocked_range<int>(0, Y.n_cols, 100), 
    [&](const tbb::blocked_range<int>& r){ 

    disable_parallel_blas();

    // iterate through responses 
    for (int j = r.begin(); j != r.end(); ++j) { 

      T1 y = Y.col(j);
      vec w = Weights.col(j);

      spectralDecomp<T3> dcmp;
      dcmp.initWithIndicator(Z, w);

      fastlmm fit = fastlmm<T1, T2, T3>(y, X, dcmp.get_vectors(), dcmp.get_values(), w, md, REML);

      fit.estimate_delta( left, right, tol );

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