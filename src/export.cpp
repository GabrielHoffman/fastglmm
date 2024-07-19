#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <omp.h>
#include "fastlmm.hpp"
#include "nb_theta.hpp"

using namespace Rcpp; 
using namespace arma;
using namespace fastglmm;

// Depends on Rcpp::List, so define outside of class
const List toList(fastlmm_result &res){

  return List::create( 
                Named("logLik") = res.logLik, 
                Named("beta")   = res.beta,
                Named("beta_se")= res.beta_se,
                Named("vcov")   = res.vcov, 
                Named("delta")  = res.delta,
                Named("sig_g")  = res.sig_g,
                Named("sig_e")  = res.sig_e,
                Named("iter")   = res.iter);
}

template <typename T1, typename T2, typename T3>
const List toList(fastlmm<T1, T2, T3> & fit){
  fastlmm_result a = fit.get_result();
  return toList( a);
}


List toList( const vector<fastlmm_result> &resList){
  List L = List::create();

  for(int i=0; i<resList.size(); i++){
    fastlmm_result a = resList.at(i);
    L.push_back( toList(a) );
  }
  return L;
}


// Cannot exprot template functions
// so write separately for matrix and sparse matrix

// [[Rcpp::export(".fastlmm_mat")]]
List fastlmm_mat( const vec &Y, 
                  const mat &X,  
                  const mat &U, 
                  const vec &s,
                  const vec &weights,
                  const double &delta,
                  const double & tol){

  // initialize
  fastlmm fit = fastlmm(Y, X, U, s, weights);

  if( delta > 0 ){
    fit.eval_delta( delta ); 
  }else{
    fit.estimate_delta( tol );
  }

  return toList(fit);
}


// [[Rcpp::export(".fastlmm_spmat")]]
List fastlmm_spmat( const vec &Y, 
                    const mat &X,  
                    const sp_mat &U, 
                    const vec &s,
                    const vec &weights,
                    const double &delta,
                    const double & tol){

  // initialize  
  fastlmm fit = fastlmm(Y, X, U, s, weights);

  if( delta > 0 ){
    fit.eval_delta( delta ); 
  }else{
    fit.estimate_delta( tol );
  }

  return toList(fit);
}


// [[Rcpp::export(".fastlmm_batch_mat")]]
List fastlmm_batch_mat( const mat &Y_all, 
                        const mat &X,  
                        const mat &U, 
                        const vec &s,
                        const mat &weights,
                        const double &delta,
                        const double & tol){

  // initialize
  fastlmm fit = fastlmm<mat, mat, mat>(X, U, s);

  vector<fastlmm_result> res;
  res = fit.fit_batch_response(Y_all, weights, delta, tol );

  return toList(res);
}


// [[Rcpp::export(".fastlmm_batch_spmat")]]
List fastlmm_batch_spmat( const mat &Y_all, 
                          const mat &X,  
                          const sp_mat &U, 
                          const vec &s,
                          const vec &weights,
                          const double &delta,
                          const double & tol){

  // initialize
  fastlmm fit = fastlmm<mat, mat, sp_mat>(X, U, s);

  vector<fastlmm_result> res;
  res = fit.fit_batch_response(Y_all, weights, delta, tol );

  return toList(res);
}



