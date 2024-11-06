#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <fastlmmLib.h>

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
  fastlmmBatchResponse fit = 
    fastlmmBatchResponse<mat, mat, mat>(Y_all, X, Z, weights, left, right, tol, nthreads);

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
  fastlmmBatchResponse fit = 
    fastlmmBatchResponse<mat, sp_mat, mat>(Y_all, X, Z, weights, left, right, tol, nthreads);

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
  fastlmmBatchResponse fit = 
    fastlmmBatchResponse<mat, mat, sp_mat>(Y_all, X, Z, weights, left, right, tol, nthreads);

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
  fastlmmBatchResponse fit = 
    fastlmmBatchResponse<mat, sp_mat, sp_mat>(Y_all, X, Z, weights, left, right, tol, nthreads);

  // evaluate each response
  vector<fastlmm_result> res = fit.eval();

  return toList(res);
}









// [[Rcpp::export(".fastlmm_batch_design_m")]]
List fastlmm_batch_design_m(const arma::mat &Y, 
                            const arma::mat &X,
                            const arma::mat &X_add,   
                            const arma::mat &U, 
                            const arma::vec &s,
                            const arma::vec &weights,
                            const double &delta,
                            const double &left,
                            const double &right,
                            const double &tol,
                            const int &nthreads){

  // initialize
  fastlmmBatchDesign fit = 
    fastlmmBatchDesign<mat, mat, mat>(Y, X, U, s, weights);

  // evaluate each column of X_add, one at a time
  vector<fastlmm_result> res = fit.eval(X_add, delta, left, right, tol, nthreads );

  return toList( res );
}


// [[Rcpp::export(".fastlmm_batch_design_s")]]
List fastlmm_batch_design_s(const arma::mat &Y, 
                            const arma::mat &X,
                            const arma::mat &X_add,   
                            const arma::sp_mat &U, 
                            const arma::vec &s,
                            const arma::vec &weights,
                            const double &delta,
                            const double &left,
                            const double &right,
                            const double &tol,
                            const int &nthreads){

  // initialize
  fastlmmBatchDesign fit = 
    fastlmmBatchDesign<mat, mat, sp_mat>(Y, X, U, s, weights);

  // evaluate each column of X_add, one at a time
  vector<fastlmm_result> res = fit.eval(X_add, delta, left, right, tol, nthreads );

  return toList( res );
}







