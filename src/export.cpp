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
#include "glmmFitFeatures.h"
#include "glmmFitResponses.h"
#include "log_moments_nb.h"

using namespace Rcpp; 
using namespace arma;
using namespace fastglmmLib;


// X = mat
// U = mat
// [[Rcpp::export(".fastlmm_mm")]]
List fastlmm_mm( 
  const arma::vec &y, 
  const arma::mat &X,  
  const arma::mat &U, 
  const arma::vec &s,
  const arma::vec &weights,
  const std::string &dcmpMethod,
  const double &delta,
  const double &left,
  const double &right,
  const double &tol,
  const double &lambda,
  const int &nthreads,
  const bool REML = false){

  ModelDetail md = MAX;

  ZTYPE type = dcmpMethod == "categorical" ? 
                CATEGORICAL : GENERAL;
  spectralDecomp dcmp(U, s, type);

  // initialize
  fastlmm fit = fastlmm(y, X, dcmp, weights, md, lambda, REML);

  if( delta > 0 ){
    fit.eval_delta( delta ); 
  }else{
    fit.estimate_delta( left, right, tol );
  }

  return toList(fit);
}

// [[Rcpp::export(".fastlmm_ms")]]
List fastlmm_ms( 
  const arma::vec &y, 
  const arma::mat &X,  
  const arma::sp_mat &U, 
  const arma::vec &s,
  const arma::vec &weights,
  const std::string &dcmpMethod,
  const double &delta,
  const double &left,
  const double &right,
  const double &tol,
  const double &lambda,
  const int &nthreads,
  const bool REML = false){

  ModelDetail md = MAX;

  ZTYPE type = dcmpMethod == "categorical" ? 
                CATEGORICAL : GENERAL;
  spectralDecomp dcmp(U, s, type);

  // initialize
  fastlmm fit = fastlmm(y, X, dcmp, weights, md, lambda, REML);

  if( delta > 0 ){
    fit.eval_delta( delta ); 
  }else{
    fit.estimate_delta( left, right, tol );
  }

  return toList(fit);
}



// [[Rcpp::export(".fastlmm_sm")]]
List fastlmm_sm( 
  const arma::vec &y, 
  const arma::sp_mat &X,  
  const arma::mat &U, 
  const arma::vec &s,
  const arma::vec &weights,
  const std::string &dcmpMethod,
  const double &delta,
  const double &left,
  const double &right,
  const double &tol,
  const double &lambda,
  const int &nthreads,
  const bool REML = false){

  ModelDetail md = MAX;

  ZTYPE type = dcmpMethod == "categorical" ? 
                CATEGORICAL : GENERAL;
  spectralDecomp dcmp(U, s, type);

  // initialize
  fastlmm fit = fastlmm(y, X, dcmp, weights, md, lambda, REML);

  if( delta > 0 ){
    fit.eval_delta( delta ); 
  }else{
    fit.estimate_delta( left, right, tol );
  }

  return toList(fit);
}

// [[Rcpp::export(".fastlmm_ss")]]
List fastlmm_ss( 
  const arma::vec &y, 
  const arma::sp_mat &X,  
  const arma::sp_mat &U, 
  const arma::vec &s,
  const arma::vec &weights,
  const std::string &dcmpMethod,
  const double &delta,
  const double &left,
  const double &right,
  const double &tol,
  const double &lambda,
  const int &nthreads,
  const bool REML = false){

  ModelDetail md = MAX;

  ZTYPE type = dcmpMethod == "categorical" ? 
                CATEGORICAL : GENERAL;
  spectralDecomp dcmp(U, s, type);

  // initialize
  fastlmm fit = fastlmm(y, X, dcmp, weights, md, lambda, REML);

  if( delta > 0 ){
    fit.eval_delta( delta ); 
  }else{
    fit.estimate_delta( left, right, tol );
  }

  return toList(fit);
}



// [[Rcpp::export(".nb_theta")]]
double nb_theta(
  const NumericVector &y,
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
List fastglmm_mm( 
  const arma::vec &y, 
  const arma::mat &X,  
  const arma::mat &U, 
  const arma::vec &s,
  const arma::vec &weights,
  const arma::vec &offset,
  const std::string &family,
  const std::string &dcmpMethod,
  const double &delta,
  const double &left,
  const double &right,
  const double &tol,
  const double &tol_eta,
  const int &maxit,
  const double &lambda,
  const int &nthreads, 
  const bool &doCoxReid = true){

  ModelDetail md = MAX;

  ZTYPE type = dcmpMethod == "categorical" ? 
                CATEGORICAL : GENERAL;
  spectralDecomp dcmp(U, s, type);

  // initialize
  fastglmm fit = fastglmm<vec,mat,mat>(y, X, dcmp, weights, offset, family, md, tol, tol_eta, maxit, lambda, delta, left, right, true, doCoxReid);

  return toList(fit);
}


// y = vec
// X = mat
// U = sp_mat
// [[Rcpp::export(".fastglmm_ms")]]
List fastglmm_ms( 
  const arma::vec &y, 
  const arma::mat &X,  
  const arma::sp_mat &U, 
  const arma::vec &s,
  const arma::vec &weights,
  const arma::vec &offset,
  const std::string &family,
  const std::string &dcmpMethod,
  const double &delta,
  const double &left,
  const double &right,
  const double &tol,
  const double &tol_eta,
  const int &maxit,
  const double &lambda,
  const int &nthreads, 
  const bool &doCoxReid = true){

  ModelDetail md = MAX;

  ZTYPE type = dcmpMethod == "categorical" ? 
                CATEGORICAL : GENERAL;
  spectralDecomp dcmp(U, s, type);

  // initialize
  fastglmm fit = fastglmm<vec,mat,sp_mat>(y, X, dcmp, weights, offset, family, md, tol, tol_eta, maxit, lambda, delta, left, right, true, doCoxReid);

  return toList(fit);
}


/* Log Moments of NB given mu vector
*/
// [[Rcpp::export]]
DataFrame log_moments_nb_mu(
  const vec &mu,
  const double & theta, // overdispersion
  const double & c = 1.0, // pseudocount
  const double & p_tail = 1e-4) {

  auto [signal, noise] = _log_moments_nb_mu(mu, theta, c, p_tail);

  return DataFrame::create(
    Named("var.signal") = signal,
    Named("var.noise")  = noise
  );
}


// /* Log Moments of NB given mu matrix
// */
// // [[Rcpp::export]]
// List log_moments_nb_mat(
//   const mat & Mu,
//   const vec & theta, // overdispersion
//   const double & c = 1.0, // pseudocount
//   const double & p_tail = 1e-4) {

//   int n_responses = Mu.n_cols;
//   vec signal(n_responses);
//   vec noise(n_responses);

//   for(int i=0; i<n_responses; ++i){
//     auto [signal_, noise_] = _log_moments_nb_mu(Mu.col(i), theta(i), c, p_tail);    
//     signal(i) = signal_;
//     noise(i) = noise_; 
//   }

//   return List::create(
//     Named("var.signal") = signal,
//     Named("var.noise")  = noise
//   );
// }


/* Log Moments of NB given X, Beta, and offset */
// [[Rcpp::export]]
DataFrame log_moments_nb_XB(
  const mat & X,
  const mat & Beta,
  const vec & offset,
  const vec & theta, // overdispersion
  const string &method,
  const double & c = 1.0, // pseudocount
  const double & p_tail = 1e-4,
  const int & nthreads = 10) {

  auto [signal, noise, alpha] = _log_moments_nb_XB(X, Beta, offset, theta, method, c, p_tail, nthreads);

  return DataFrame::create(
    Named("var.signal") = wrap(signal),
    Named("var.noise")  = wrap(noise),
    Named("alpha")  = wrap(alpha)
  );
}

// [[Rcpp::export]]
DataFrame log_moments_nb_BlupXBZ(
  const mat & BLUP,
  const mat & X,
  const mat & Beta,
  const mat & Z,  
  const vec & weights,
  const vec & offset,
  const vec & theta, // overdispersion
  const vec & delta, // variance component ratio
  const string &method,
  const string &dcmpMethod, 
  const double & c = 1.0, // pseudocount
  const double & p_tail = 1e-4,
  const int & nthreads = 10) {

  auto [signal, noise, alpha] = _log_moments_nb_BlupXBZ(BLUP, X, Beta, Z, weights, offset, theta, delta, method, dcmpMethod, c, p_tail, nthreads);

  return DataFrame::create(
    Named("var.signal") = wrap(signal),
    Named("var.noise")  = wrap(noise),
    Named("alpha")  = wrap(alpha)
  );
}

// [[Rcpp::export]]
DataFrame log_moments_nb_BlupXBZ_sp(
  const mat & BLUP,
  const mat & X,
  const mat & Beta,
  const sp_mat & Z,  
  const vec & weights,
  const vec & offset,
  const vec & theta, // overdispersion
  const vec & delta, // variance component ratio
  const string &method,
  const string &dcmpMethod, 
  const double & c = 1.0, // pseudocount
  const double & p_tail = 1e-4,
  const int & nthreads = 10) {

   auto [signal, noise, alpha] = _log_moments_nb_BlupXBZ(BLUP, X, Beta, Z, weights, offset, theta, delta, method, dcmpMethod, c, p_tail, nthreads);

  return DataFrame::create(
    Named("var.signal") = wrap(signal),
    Named("var.noise")  = wrap(noise),
    Named("alpha")  = wrap(alpha)
  );
}




