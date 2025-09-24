#ifndef GLMM_FIT_FEATURES_H_
#define GLMM_FIT_FEATURES_H_

#include "fastglmm_fit.h"
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
class glmmFitFeatures {

  public:
  glmmFitFeatures( const T1 &y, 
                  const T2 &X, 
                  const string &family, 
                  const T3 &U,
                  const vec &s,
                  const vec &weights = {},
                  const vec &offset = {},
                  const double &delta = -1,
                  const double &left = -10,
                  const double &right = 10,
                  const double &tol = 1e-5,
                  const double &tol_eta = 1e-5,
                  const int &nthreads = 1,
                  const ModelDetail md = LOW);

  ModelFitGLMMList eval(const T2 &X_add_,
                      const vector<string> &ids);

  private:
  T1 y; 
  T2 X_shared;  
  string family;
  T3 U;
  vec s, weights, offset;
  double delta, left, right, tol, tol_eta;
  int nthreads;
  ModelDetail md;
  spectralDecomp<T3> dcmp;
};



// constructor
template <typename T1, typename T2, typename T3> 
glmmFitFeatures<T1, T2, T3>::glmmFitFeatures(
                      const T1 &y, 
                      const T2 &X, 
                      const string &family, 
                      const T3 &U,
                      const vec &s,
                      const vec &weights, 
                      const vec &offset,    
                      const double &delta,
                      const double &left,
                      const double &right,
                      const double &tol,
                      const double &tol_eta,
                      const int &nthreads,
                      const ModelDetail md) :
  y(y), X_shared(X), family(family), weights(weights), offset(offset),
  delta(delta), left(left), right(right), tol(tol), tol_eta(tol_eta), nthreads(nthreads), md(md)
  {

  dcmp.initWithEigenDecomp(U, s, weights);
}



template <typename T1, typename T2, typename T3> 
ModelFitGLMMList 
glmmFitFeatures<T1, T2, T3>::eval( const T2 &X_add_,
                                  const vector<string> &ids){

  int n_tests = X_add_.n_cols;

  // store results
  ModelFitGLMMList result(n_tests, ModelFitGLMM());

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
      fastglmm fit = fastglmm(y, X_combined, dcmp.get_vectors(), dcmp.get_values(), weights, offset, family, md, tol, tol_eta);

      result.at(j) = fit.get_result();
      result.at(j).ID = ids[j];
    }
  }); });

  return result;
}



} // end namespace



#endif