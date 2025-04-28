#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "fastlmmLib.h"
#include "exportToR_fastlmm.h"

using namespace arma;
using namespace fastlmmLib;

// Cannot export template functions
// so write separately for matrix and sparse matrix

// Y = vec
// X = mat
// U = mat
// [[Rcpp::export(".fastlmm_vmm")]]
List fastlmm_vmm( const arma::vec &Y, 
                  const arma::mat &X,  
                  const arma::mat &U, 
                  const arma::vec &s,
                  const arma::vec &weights,
                  const double &delta,
                  const double &left,
                  const double &right,
                  const double &tol,
                  const int &nthreads){

  // initialize
  fastlmm fit = fastlmm(Y, X, U, s, weights);

  if( delta > 0 ){
    fit.eval_delta( delta ); 
  }else{
    fit.estimate_delta( left, right, tol );
  }

  return toList(fit);
}

// [[Rcpp::export(".fastlmm_vms")]]
List fastlmm_vms( const arma::vec &Y, 
                  const arma::mat &X,  
                  const arma::sp_mat &U, 
                  const arma::vec &s,
                  const arma::vec &weights,
                  const double &delta,
                  const double &left,
                  const double &right,
                  const double &tol,
                  const int &nthreads){

  // initialize
  fastlmm fit = fastlmm(Y, X, U, s, weights);

  if( delta > 0 ){
    fit.eval_delta( delta ); 
  }else{
    fit.estimate_delta( left, right, tol );
  }

  return toList(fit);
}



// [[Rcpp::export(".fastlmm_vsm")]]
List fastlmm_vsm( const arma::vec &Y, 
                  const arma::sp_mat &X,  
                  const arma::mat &U, 
                  const arma::vec &s,
                  const arma::vec &weights,
                  const double &delta,
                  const double &left,
                  const double &right,
                  const double &tol,
                  const int &nthreads){

  // initialize
  fastlmm fit = fastlmm(Y, X, U, s, weights);

  if( delta > 0 ){
    fit.eval_delta( delta ); 
  }else{
    fit.estimate_delta( left, right, tol );
  }

  return toList(fit);
}

// [[Rcpp::export(".fastlmm_vss")]]
List fastlmm_vss( const arma::vec &Y, 
                  const arma::sp_mat &X,  
                  const arma::sp_mat &U, 
                  const arma::vec &s,
                  const arma::vec &weights,
                  const double &delta,
                  const double &left,
                  const double &right,
                  const double &tol,
                  const int &nthreads){

  // initialize
  fastlmm fit = fastlmm(Y, X, U, s, weights);

  if( delta > 0 ){
    fit.eval_delta( delta ); 
  }else{
    fit.estimate_delta( left, right, tol );
  }

  return toList(fit);
}




// Y = mat
// X = mat
// U = mat
// [[Rcpp::export(".fastlmm_mmm")]]
List fastlmm_mmm(   const arma::mat &Y_all, 
                    const arma::mat &X,  
                    const arma::mat &Z,
                    const arma::mat &weights,
                    const double &delta,
                    const double &left,
                    const double &right,
                    const double &tol,
                    const int &nthreads){

  // initialize
  lmmFitResponses fit = 
    lmmFitResponses<mat, mat, mat>(Y_all, X, Z, weights, left, right, tol, nthreads);

  // evaluate each response
  vector<ModelFitLMM> res = fit.eval();

  return toList(res);
}


// [[Rcpp::export(".fastlmm_msm")]]
List fastlmm_msm(   const arma::mat &Y_all, 
                    const arma::sp_mat &X,  
                    const arma::mat &Z,
                    const arma::mat &weights,
                    const double &delta,
                    const double &left,
                    const double &right,
                    const double &tol,
                    const int &nthreads){

  // initialize
  lmmFitResponses fit = 
    lmmFitResponses<mat, sp_mat, mat>(Y_all, X, Z, weights, left, right, tol, nthreads);

  // evaluate each response
  vector<ModelFitLMM> res = fit.eval();

  return toList(res);
}

// [[Rcpp::export(".fastlmm_mms")]]
List fastlmm_mms(   const arma::mat &Y_all, 
                    const arma::mat &X,  
                    const arma::sp_mat &Z,
                    const arma::mat &weights,
                    const double &left,
                    const double &right,
                    const double &tol,
                    const int &nthreads){

  // initialize
  lmmFitResponses fit = 
    lmmFitResponses<mat, mat, sp_mat>(Y_all, X, Z, weights, left, right, tol, nthreads);

  // evaluate each response
  vector<ModelFitLMM> res = fit.eval();

  return toList(res);
}


// [[Rcpp::export(".fastlmm_mss")]]
List fastlmm_mss(   const arma::mat &Y_all, 
                    const arma::sp_mat &X,  
                    const arma::sp_mat &Z,
                    const arma::mat &weights,
                    const double &delta,
                    const double &left,
                    const double &right,
                    const double &tol,
                    const int &nthreads){

  // initialize
  lmmFitResponses fit = 
    lmmFitResponses<mat, sp_mat, sp_mat>(Y_all, X, Z, weights, left, right, tol, nthreads);

  // evaluate each response
  vector<ModelFitLMM> res = fit.eval();

  return toList(res);
}


