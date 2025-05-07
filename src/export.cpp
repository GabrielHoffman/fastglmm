/***************************************************************
 * @file    export.cpp
 * @author  Gabriel Hoffman
 * @email   gabriel.hoffman@mssm.edu
 * @brief   Export functions to R
 * Copyright (C) 2024 Gabriel Hoffman
 **************************************************************/

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "fastglmm.h"
#include "exportToR_fastlmm.h"
#include "nb_theta.h"

using namespace Rcpp; 
using namespace arma;
using namespace fastglmmLib;



// Cannot export template functions
// so write separately for matrix and sparse matrix

// Y = vec
// X = mat
// U = mat
// [[Rcpp::export(".fastlmm_vmm")]]
List fastlmm_vmm( const arma::vec &y, 
                  const arma::mat &X,  
                  const arma::mat &U, 
                  const arma::vec &s,
                  const arma::vec &weights,
                  const double &delta,
                  const double &left,
                  const double &right,
                  const double &tol,
                  const int &nthreads,
                  const bool REML = false){

  ModelDetail md = MAX;

  // initialize
  fastlmm fit = fastlmm(y, X, U, s, weights, md, REML);

  if( delta > 0 ){
    fit.eval_delta( delta ); 
  }else{
    fit.estimate_delta( left, right, tol );
  }

  return toList(fit);
}

// [[Rcpp::export(".fastlmm_vms")]]
List fastlmm_vms( const arma::vec &y, 
                  const arma::mat &X,  
                  const arma::sp_mat &U, 
                  const arma::vec &s,
                  const arma::vec &weights,
                  const double &delta,
                  const double &left,
                  const double &right,
                  const double &tol,
                  const int &nthreads,
                  const bool REML = false){

  ModelDetail md = MAX;

  // initialize
  fastlmm fit = fastlmm(y, X, U, s, weights, md, REML);

  if( delta > 0 ){
    fit.eval_delta( delta ); 
  }else{
    fit.estimate_delta( left, right, tol );
  }

  return toList(fit);
}



// [[Rcpp::export(".fastlmm_vsm")]]
List fastlmm_vsm( const arma::vec &y, 
                  const arma::sp_mat &X,  
                  const arma::mat &U, 
                  const arma::vec &s,
                  const arma::vec &weights,
                  const double &delta,
                  const double &left,
                  const double &right,
                  const double &tol,
                  const int &nthreads,
                  const bool REML = false){

  ModelDetail md = MAX;

  // initialize
  fastlmm fit = fastlmm(y, X, U, s, weights, md, REML);

  if( delta > 0 ){
    fit.eval_delta( delta ); 
  }else{
    fit.estimate_delta( left, right, tol );
  }

  return toList(fit);
}

// [[Rcpp::export(".fastlmm_vss")]]
List fastlmm_vss( const arma::vec &y, 
                  const arma::sp_mat &X,  
                  const arma::sp_mat &U, 
                  const arma::vec &s,
                  const arma::vec &weights,
                  const double &delta,
                  const double &left,
                  const double &right,
                  const double &tol,
                  const int &nthreads,
                  const bool REML = false){

  ModelDetail md = MAX;

  // initialize
  fastlmm fit = fastlmm(y, X, U, s, weights, md, REML);

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
List fastlmm_mmm(   const arma::mat &Y, 
                    const std::vector<std::string> ids, 
                    const arma::mat &X,  
                    const arma::mat &Z,
                    const arma::mat &Weights,
                    const double &delta,
                    const double &left,
                    const double &right,
                    const double &tol,
                    const int &nthreads,
                    const bool REML = false){

  ModelDetail md = MAX;

  // initialize
  lmmFitResponses fit = 
    lmmFitResponses<mat, mat, mat>(X, Z, left, right, tol, nthreads, md, REML);

  // evaluate each response
  vector<ModelFitLMM> res = fit.eval(Y, ids, Weights);

  return toList(res);
}


// [[Rcpp::export(".fastlmm_msm")]]
List fastlmm_msm(   const arma::mat &Y, 
                    const std::vector<std::string> ids, 
                    const arma::sp_mat &X,  
                    const arma::mat &Z,
                    const arma::mat &Weights,
                    const double &delta,
                    const double &left,
                    const double &right,
                    const double &tol,
                    const int &nthreads,
                    const bool REML = false){

  ModelDetail md = MAX;

  // initialize
  lmmFitResponses fit = 
    lmmFitResponses<mat, sp_mat, mat>(X, Z, left, right, tol, nthreads, md, REML);

  // evaluate each response
  vector<ModelFitLMM> res = fit.eval(Y, ids, Weights);

  return toList(res);
}

// [[Rcpp::export(".fastlmm_mms")]]
List fastlmm_mms(   const arma::mat &Y, 
                    const std::vector<std::string> ids, 
                    const arma::mat &X,  
                    const arma::sp_mat &Z,
                    const arma::mat &Weights,
                    const double &left,
                    const double &right,
                    const double &tol,
                    const int &nthreads,
                    const bool REML = false){

  ModelDetail md = MAX;

  // initialize
  lmmFitResponses fit = 
    lmmFitResponses<mat, mat, sp_mat>(X, Z, left, right, tol, nthreads, md, REML);

  // evaluate each response
  vector<ModelFitLMM> res = fit.eval(Y, ids, Weights);

  return toList(res);
}


// [[Rcpp::export(".fastlmm_mss")]]
List fastlmm_mss(   const arma::mat &Y,
                    const std::vector<std::string> ids, 
                    const arma::sp_mat &X,  
                    const arma::sp_mat &Z,
                    const arma::mat &Weights,
                    const double &delta,
                    const double &left,
                    const double &right,
                    const double &tol,
                    const int &nthreads,
                    const bool REML = false){

  ModelDetail md = MAX;

  // initialize
  lmmFitResponses fit = 
    lmmFitResponses<mat, sp_mat, sp_mat>(X, Z, left, right, tol, nthreads, md, REML);

  // evaluate each response
  vector<ModelFitLMM> res = fit.eval(Y, ids, Weights);

  return toList(res);
}


// [[Rcpp::export(".nb_theta")]]
double nb_theta(const NumericVector &y,
                const NumericVector &mu,
                const double &n, 
                const NumericVector &weights,
                const double &left = -5,
                const double &right = 20,
                const double &tol = 1e-5){

  return nb_theta_ml(y, mu, n, weights, {}, false, {}, left, right, tol);
}


// y = vec
// X = mat
// U = mat
// [[Rcpp::export(".fastglmm_mm")]]
List fastglmm_mm( const arma::vec &y, 
                  const arma::mat &X,  
                  const arma::mat &U, 
                  const arma::vec &s,
                  const arma::vec &weights,
                  const arma::vec &offset,
                  const std::string &family,
                  const double &delta,
                  const double &left,
                  const double &right,
                  const double &tol,
                  const int &nthreads){

  ModelDetail md = MAX;
  double tol_eta = 1e-4;

  // initialize
  fastglmm fit = fastglmm<vec,mat,mat>(y, X, U, s, weights, offset, family, md, tol, tol_eta, true);

  return toList(fit);
}


// y = vec
// X = mat
// U = mat
// [[Rcpp::export(".fastglmm_ms")]]
List fastglmm_ms( const arma::vec &y, 
                  const arma::mat &X,  
                  const arma::sp_mat &U, 
                  const arma::vec &s,
                  const arma::vec &weights,
                  const arma::vec &offset,
                  const std::string &family,
                  const double &delta,
                  const double &left,
                  const double &right,
                  const double &tol,
                  const int &nthreads){

  ModelDetail md = MAX;
  double tol_eta = 1e-4;

  // initialize
  fastglmm fit = fastglmm<vec,mat,sp_mat>(y, X, U, s, weights, offset, family, md, tol, tol_eta, true);

  return toList(fit);
}


// [[Rcpp::export]]
void test_time_m(const arma::mat &U, 
                const arma::vec &s,
                const arma::vec &weights){

  spectralDecomp<mat> dcmp;
  mat Z = scaleEachRow(U, sqrt(s));
  vec s2;
  mat U2;

  // for(int i=0; i<300; i++){    
  //   // s2 = (weights.t() * Z).t();
  //   mat Z = scaleEachRow(U, sqrt(s));
  // }

  for(int i=0; i<300; i++){    
    U2 = scaleRowsCols(Z, (weights), s);
  }
}



// [[Rcpp::export]]
void test_time_s(const arma::sp_mat &U, 
                const arma::vec &s,
                const arma::vec &weights){

  spectralDecomp<sp_mat> dcmp;
  sp_mat Z = scaleEachRow(U, sqrt(s));
  vec s2;
  sp_mat U2;  
  
  for(int i=0; i<300; i++){    
    U2 = scaleRowsCols(Z, (weights), s);
  }

  //   s2 = (weights.t() * Z).t();
  // // sp_mat a = sp_mat(arma::diagmat(weights));
  // sp_mat a(weights.n_elem, weights.n_elem);
  // a.diag() = weights;
  // a.brief_print();
  // // sp_mat b = sp_mat(arma::diagmat(s2));
  // sp_mat b(s2.n_elem, s2.n_elem);
  // b.diag() = s2;

  // for(int i=0; i<300; i++){    
  //   U2 = a * Z * b;
  // }
}








