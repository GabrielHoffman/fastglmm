#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <omp.h>
#include "fastlmm.hpp"
#include "nb_theta.hpp"

using namespace Rcpp; 
using namespace arma;

// Depends on Rcpp::List, so define outside of class
const List toList(FASTLMM_result &res){

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

template <typename T>
const List toList(FASTLMM<T> & fit){
  FASTLMM_result a = fit.get_result();
  return toList( a);
}


List toList( const std::vector<FASTLMM_result> &resList){
  List L = List::create();

  for(int i=0; i<resList.size(); i++){
    FASTLMM_result a = resList.at(i);
    L.push_back( toList(a) );
  }
  return L;
}


// Cannot exprot template functions
// so write separately for matrix and sparse matrix

// [[Rcpp::export(".fastlmm_mat")]]
List fastlmm_mat( const arma::vec &Y, 
                  const arma::mat &X,  
                  const arma::mat &U, 
                  const arma::vec &s,
                  const arma::vec &weights,
                  const double &delta,
                  const double & tol){

  // initialize
  FASTLMM<arma::mat> fit = 
      FASTLMM<arma::mat>(Y, X, U, s, weights);

  if( delta > 0 ){
    fit.eval_delta( delta ); 
  }else{
    fit.estimate_delta( tol );
  }

  return toList(fit);
}


// [[Rcpp::export(".fastlmm_spmat")]]
List fastlmm_spmat( const arma::vec &Y, 
                    const arma::mat &X,  
                    const arma::sp_mat &U, 
                    const arma::vec &s,
                    const arma::vec &weights,
                    const double &delta,
                    const double & tol){

  // initialize  
  FASTLMM<arma::sp_mat> fit = 
      FASTLMM<arma::sp_mat>(Y, X, U, s, weights);

  if( delta > 0 ){
    fit.eval_delta( delta ); 
  }else{
    fit.estimate_delta( tol );
  }

  return toList(fit);
}


// [[Rcpp::export(".fastlmm_batch_mat")]]
List fastlmm_batch_mat( const arma::mat &Y_all, 
                        const arma::mat &X,  
                        const arma::mat &U, 
                        const arma::vec &s,
                        const arma::mat &weights,
                        const double &delta,
                        const double & tol){
  // initialize
  FASTLMM<arma::mat> fit = 
      FASTLMM<arma::mat>(X, U, s);

  std::vector<FASTLMM_result> res;
  res = fit.fit_batch_response(Y_all, weights, delta, tol );

  return toList(res);
}


// [[Rcpp::export(".fastlmm_batch_spmat")]]
List fastlmm_batch_spmat( const arma::mat &Y_all, 
                          const arma::mat &X,  
                          const arma::sp_mat &U, 
                          const arma::vec &s,
                          const arma::vec &weights,
                          const double &delta,
                          const double & tol){
  // initialize
  FASTLMM<arma::sp_mat> fit = 
      FASTLMM<arma::sp_mat>(X, U, s);

  std::vector<FASTLMM_result> res;
  res = fit.fit_batch_response(Y_all, weights, delta, tol );

  return toList(res);
}






