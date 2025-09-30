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

using namespace Rcpp; 
using namespace arma;
using namespace fastglmmLib;


// X = mat
// U = mat
// [[Rcpp::export(".fastlmm_mm")]]
List fastlmm_mm( const arma::vec &y, 
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

  spectralDecomp dcmp(U, s);

  // initialize
  fastlmm fit = fastlmm(y, X, dcmp, weights, md, REML);

  if( delta > 0 ){
    fit.eval_delta( delta ); 
  }else{
    fit.estimate_delta( left, right, tol );
  }

  return toList(fit);
}

// [[Rcpp::export(".fastlmm_ms")]]
List fastlmm_ms( const arma::vec &y, 
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

  spectralDecompCategorical dcmp(U, s);

  // initialize
  fastlmm fit = fastlmm(y, X, dcmp, weights, md, REML);

  if( delta > 0 ){
    fit.eval_delta( delta ); 
  }else{
    fit.estimate_delta( left, right, tol );
  }

  return toList(fit);
}



// [[Rcpp::export(".fastlmm_sm")]]
List fastlmm_sm( const arma::vec &y, 
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

  spectralDecomp dcmp(U, s);

  // initialize
  fastlmm fit = fastlmm(y, X, dcmp, weights, md, REML);

  if( delta > 0 ){
    fit.eval_delta( delta ); 
  }else{
    fit.estimate_delta( left, right, tol );
  }

  return toList(fit);
}

// [[Rcpp::export(".fastlmm_ss")]]
List fastlmm_ss( const arma::vec &y, 
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

  spectralDecompCategorical dcmp(U, s);

  // initialize
  fastlmm fit = fastlmm(y, X, dcmp, weights, md, REML);

  if( delta > 0 ){
    fit.eval_delta( delta ); 
  }else{
    fit.estimate_delta( left, right, tol );
  }

  return toList(fit);
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
                  const double &tol_eta,
                  const int &maxit,
                  const int &nthreads){

  ModelDetail md = MAX;

  spectralDecomp dcmp(U, s);

  // initialize
  fastglmm fit = fastglmm<vec,mat,mat>(y, X, dcmp, weights, offset, family, md, tol, tol_eta, maxit, delta, left, right, true);

  return toList(fit);
}


// y = vec
// X = mat
// U = sp_mat
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
                  const double &tol_eta,
                  const int &maxit,
                  const int &nthreads){

  ModelDetail md = MAX;

  spectralDecompCategorical dcmp(U, s);

  // initialize
  fastglmm fit = fastglmm<vec,mat,sp_mat>(y, X, dcmp, weights, offset, family, md, tol, tol_eta, maxit, delta, left, right, true);

  return toList(fit);
}


