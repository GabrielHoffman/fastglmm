#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "fastlmmLib.h"
#include "linearRegression.h"
#include "exportToR.h"

using namespace Rcpp; 
using namespace arma;
using namespace fastlmmLib;




// Depends on Rcpp::List, so define outside of class
const List toList(fastlmm_result &res){

  return List::create( 
                Named("logLik")       = res.logLik, 
                Named("coefficients") = res.beta,
                Named("se")           = res.beta_se,
                Named("vcov")         = res.vcov, 
                Named("weights")      = res.weights,
                Named("delta")        = res.delta,
                Named("sigSq_g")      = res.sigSq_g,
                Named("sigSq_e")      = res.sigSq_e,
                Named("ru")           = res.ru,
                Named("y")            = res.y,
                Named("iter")         = res.iter);
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
  fastlmmFitResponses fit = 
    fastlmmFitResponses<mat, mat, mat>(Y_all, X, Z, weights, left, right, tol, nthreads);

  // evaluate each response
  vector<fastlmm_result> res = fit.eval();

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
  fastlmmFitResponses fit = 
    fastlmmFitResponses<mat, sp_mat, mat>(Y_all, X, Z, weights, left, right, tol, nthreads);

  // evaluate each response
  vector<fastlmm_result> res = fit.eval();

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
  fastlmmFitResponses fit = 
    fastlmmFitResponses<mat, mat, sp_mat>(Y_all, X, Z, weights, left, right, tol, nthreads);

  // evaluate each response
  vector<fastlmm_result> res = fit.eval();

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
  fastlmmFitResponses fit = 
    fastlmmFitResponses<mat, sp_mat, sp_mat>(Y_all, X, Z, weights, left, right, tol, nthreads);

  // evaluate each response
  vector<fastlmm_result> res = fit.eval();

  return toList(res);
}


// [[Rcpp::export(".fastlmmFitFeatures_m")]]
List fastlmmFitFeatures_m(const arma::mat &Y, 
                            const arma::mat &X_design,
                            const arma::mat &X_features,   
                            const arma::mat &U, 
                            const arma::vec &s,
                            const arma::vec &weights,
                            const double &delta,
                            const double &left,
                            const double &right,
                            const double &tol,
                            const int &nthreads){

  // initialize
  fastlmmFitFeatures fit = 
    fastlmmFitFeatures<mat, mat, mat>(Y, X_design, U, s, weights);

  // evaluate each column of X_add, one at a time
  vector<fastlmm_result> res = fit.eval(X_features, delta, left, right, tol, nthreads );

  return toList( res );
}


// [[Rcpp::export(".fastlmm_batch_design_s")]]
List fastlmmFitFeatures_s(const arma::mat &Y, 
                            const arma::mat &X_design, 
                            const arma::mat &X_features,  
                            const arma::sp_mat &U, 
                            const arma::vec &s,
                            const arma::vec &weights,
                            const double &delta,
                            const double &left,
                            const double &right,
                            const double &tol,
                            const int &nthreads){

  // initialize
  fastlmmFitFeatures fit = 
    fastlmmFitFeatures<mat, mat, sp_mat>(Y, X_design, U, s, weights);

  // evaluate each column of X_add, one at a time
  vector<fastlmm_result> res = fit.eval(X_features, delta, left, right, tol, nthreads );

  return toList( res );
}








// [[Rcpp::export]]
List lmFitFeatures_export(const arma::vec &y, 
                          const arma::mat &X_design, 
                          const arma::mat &X_features, 
                          const vector<string> &ids, 
                          const arma::vec &weights, 
                          const int detail = 0, 
                          const bool &preprojection = true, 
                          const int &nthreads = 1){

  ModelDetail md = static_cast<ModelDetail>(detail);

  if( preprojection && md == MOST ){
    stop("Cannot compute hatvalues with pre-projection");
  }

  vector<ModelFit> fitList = lmFitFeatures(y, X_design, X_features, ids, weights, md, preprojection, nthreads);

  return toList(fitList);
}

// [[Rcpp::export]]
List lmFitResponses_export(const arma::mat &Y, const arma::mat &X, const vector<string> &ids, const arma::mat &Weights, const int detail = 0,const int &nthreads = 1){

  ModelDetail md = static_cast<ModelDetail>(detail);


  // convert responses from __rows__ to __columns__
  vector<ModelFit> fitList = lmFitResponses(Y.t(), X, ids, Weights.t(), md, nthreads);

  return toList(fitList);
}

