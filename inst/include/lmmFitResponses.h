#ifndef LMM_FIT_RESPONSE_H_
#define LMM_FIT_RESPONSE_H_

// if -D USE_R, use RcppArmadillo library
#ifdef USE_R
// [[Rcpp::depends(RcppParallel)]]  
#include <RcppArmadillo.h>
#else
#include <armadillo>
#endif

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
                  const spectralDecomp<T3> &dcmp,
                  const double &left = -10,
                  const double &right = 10,
                  const double &tol = 1e-6,
                  const double &lambda = 0,
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
  spectralDecomp<T3> dcmp;
  double left, right, tol, lambda;
  int nthreads;
  ModelDetail md;
  bool REML;
  uvec idx_drop;
  T2 X_clean;
};



// constructor
template <typename T1, typename T2, typename T3> 
lmmFitResponses<T1, T2, T3>::lmmFitResponses(
        const T2 &X, 
        const spectralDecomp<T3> &dcmp,
        const double &left,
        const double &right,
        const double &tol,
        const double &lambda,
        const int &nthreads,
        const ModelDetail md,
        const bool REML):
  X(X), 
  dcmp(dcmp),
  left(left),
  right(right),
  tol(tol),
  lambda(lambda),
  nthreads(nthreads),
  md(md),
  REML(REML) {

  // find rows in X with NAN values
  idx_drop = rows_with_nan( X );  
  X_clean = X;
  rows_to_zero(X_clean, idx_drop);
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
    tbb::blocked_range<int>(0, Y.n_cols, 10), 
    [&](const tbb::blocked_range<int>& r){ 

    disable_parallel_blas();

    // local workspace for thread
    T1 y;
    vec w;
    uvec idx;

    // iterate through responses 
    for (int j = r.begin(); j != r.end(); ++j) { 

      y = Y.col(j);
      w = Weights.col(j);

      idx = unique(join_cols(find_nan(y), idx_drop));
      y.elem(idx).zeros();
      w.elem(idx).zeros();

      fastlmm fit = fastlmm<T1, T2, T3>(y, X_clean, dcmp, w, md, lambda, REML);

      fit.estimate_delta( left, right, tol );

      result.at(j) = fit.get_result();
      result.at(j).ID = ids[j];
    }
  }); });

  return result;
}


}




#endif