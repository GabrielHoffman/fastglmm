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
    glmmFitResponses(const T2 &X, 
                    const T3 &U,
                    const vec &s,
                    const vec &weights = {}, 
                    const vec &offset = {}, 
                    const double &left = -10,
                    const double &right = 10,
                    const double &tol = 1e-5,
                    const double &tol_eta = 1e-5,
                    const int &nthreads = 1,
                    const ModelDetail md = LOW);

    ModelFitGLMMList eval(
                    const T1 &Y,
                    const vector<string> &ids,
                    const vector<string> &family);

    private:
    T2 X;  
    T3 U;
    vec s, weights, offset;
    double left, right, tol, tol_eta;
    int nthreads;
    ModelDetail md;
};



// constructor
template <typename T1, typename T2, typename T3> 
glmmFitResponses<T1, T2, T3>::glmmFitResponses(
                            const T2 &X, 
                            const T3 &U,
                            const vec &s,
                            const vec &weights, 
                            const vec &offset, 
                            const double &left,
                            const double &right,
                            const double &tol,
                            const double &tol_eta,
                            const int &nthreads,
                            const ModelDetail md):
            X(X), 
            U(U),
            s(s),
            weights(weights),
            offset(offset),
            left(left),
            right(right),
            tol(tol),
            tol_eta(tol_eta),
            nthreads(nthreads),
            md(md) {
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

        // iterate through responses 
        for (int j = r.begin(); j != r.end(); ++j) { 

            vec y = Y.col(j);

            fastglmm fit = fastglmm<vec, T2, T3>(y, X, U, s, weights, offset, family[j], md, tol, tol_eta);

            result.at(j) = fit.get_result();
            result.at(j).ID = ids[j];
        }
    }); });
  
    return result;
}


} // end namespace


#endif