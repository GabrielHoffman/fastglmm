#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "fastlmm.hpp"

using namespace Rcpp; 
using namespace arma;

// [[Rcpp::export]]
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

  return List::create(Named("logLik") = fit.get_logLik(), 
                      Named("beta")   = fit.get_beta(),
                      Named("beta_se")   = fit.get_beta_se(),
                      Named("vcov")   = fit.get_vcov(), 
                      Named("delta")   = fit.get_delta(),
                      Named("sigg")   = fit.get_sigg(),
                      Named("sige")   = fit.get_sige(),
                      Named("iter")  = fit.get_iter());
}


// // [[Rcpp::export]]
// List fastlmm_batchY_c( const arma::mat &Y_all, 
//                   const arma::mat &X,  
//                   const arma::mat &U, 
//                   const arma::vec &s){

//   arma::mat Yu_all = U * Y_all;

//   // initialize
//   FASTLMM fit = FASTLMM(Y_all.row(0), X, U, s);

//   // iterate thru responses i.e. rows
//   for( int i = 0; i < Y_all.n_rows; i++){
//     // fit.update_Y_Yu(Y_all.row(i), Yu_all.row(i));

//     fit.eval_delta(1.0);

//     // save results

//   }

//   return List::create(Named("logLik") = fit.get_logLik(), 
//                Named("beta") = fit.get_beta());

// }


// List fastlmm_batchX_c()