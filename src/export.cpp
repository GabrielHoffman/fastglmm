#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <omp.h>
#include "fastlmm.hpp"

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

const List toList(FASTLMM & fit){
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


// [[Rcpp::export(".fastlmm_c")]]
List fastlmm_c( const arma::vec &Y, 
                  const arma::mat &X,  
                  const arma::mat &U, 
                  const arma::vec &s,
                  const double &delta){
  // initialize
  FASTLMM fit = FASTLMM(Y, X, U, s);

  if( delta > 0 ){
    fit.eval_delta( delta ); 
  }else{
    fit.estimate_delta();
  }

  return toList(fit);
}


// [[Rcpp::export(".fastlmm_batch_c")]]
List fastlmm_batch_c( const arma::mat &Y_all, 
                  const arma::mat &X,  
                  const arma::mat &U, 
                  const arma::vec &s,
                  const double &delta){
  // initialize
  FASTLMM fit = FASTLMM(X, U, s);

  std::vector<FASTLMM_result> res;
  res = fit.fit_batch_response(Y_all, delta);

  return toList(res);
}



  // List a = toList(result);