#ifndef GLMM_FIT_RESPONSE_H_
#define GLMM_FIT_RESPONSE_H_

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

#include "fastglmm_fit.h"
#include "ModelFit.h"

using namespace arma;
using namespace std;

namespace fastglmmLib {

// Order of template variables
// T1 Y
// T2 X
// T3 Z
template <typename T1, typename T2, typename T3> 
class glmmFitResponses {

  public:
  glmmFitResponses(
    const T2 &X, 
    const spectralDecomp<T3> &dcmp,
    const vec &weights = {}, 
    const vec &offset = {}, 
    const double &left = -10,
    const double &right = 10,
    const double &tol = 1e-5,
    const double &tol_eta = 1e-7,
    const int &maxit = 100,
    const double &lambda = 0,
    const int &nthreads = 1,
    const ModelDetail md = LOW);

  ModelFitGLMMList eval(
    const T1 &Y,
    const vector<string> &ids,
    const vector<string> &family);

  private:
  T2 X;  
  spectralDecomp<T3> dcmp;
  vec weights, offset;
  int maxit;
  double left, right, tol, tol_eta, lambda;
  int nthreads;
  ModelDetail md;

  uvec idx_drop;
  T2 X_clean;
};



// constructor
template <typename T1, typename T2, typename T3> 
glmmFitResponses<T1, T2, T3>::glmmFitResponses(
      const T2 &X, 
      const spectralDecomp<T3> &dcmp,
      const vec &weights, 
      const vec &offset, 
      const double &left,
      const double &right,
      const double &tol,
      const double &tol_eta,
      const int &maxit,
      const double &lambda,
      const int &nthreads,
      const ModelDetail md):
  X(X), 
  dcmp(dcmp),
  weights(weights),
  offset(offset),
  left(left),
  right(right),
  tol(tol),
  tol_eta(tol_eta),
  maxit(maxit),
  lambda(lambda),
  nthreads(nthreads),
  md(md) {

  // find rows in X with NAN values
  idx_drop = rows_with_nan(X);  
  X_clean = X;
  rows_to_zero(X_clean, idx_drop);
}



template <typename T1, typename T2, typename T3> 
ModelFitGLMMList 
  glmmFitResponses<T1, T2, T3>::eval(
      const T1 &Y,
      const vector<string> &ids,
      const vector<string> &family){

  // store results
  ModelFitGLMMList result(Y.n_cols, ModelFitGLMM());

  // Parallel part using Thread Building Blocks
  tbb::task_arena limited_arena(nthreads);
  limited_arena.execute([&] {
  tbb::parallel_for(
    tbb::blocked_range<int>(0, Y.n_cols, 100), 
    [&](const tbb::blocked_range<int>& r){ 

    disable_parallel_blas();

    T1 y;
    vec w;
    uvec idx;

    // iterate through responses 
    for (int j = r.begin(); j != r.end(); ++j) { 

      y = Y.col(j);
      w = weights;

      idx = unique(join_cols(find_nan(y), idx_drop));
      y.elem(idx).zeros();
      w.elem(idx).zeros();

      fastglmm fit = fastglmm<vec, T2, T3>(y, X_clean, dcmp, w, offset, family[j], md, tol, tol_eta, maxit, lambda);

      result.at(j) = fit.get_result();
      result.at(j).ID = ids[j];
    }
  }); });
  
  return result;
}


} // end namespace


#endif