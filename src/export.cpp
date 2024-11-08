#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <fastlmmLib.h>
#include "linearRegression.h"

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



List toList( const vector<ModelFit> & fitList){

    int ncoef = fitList[0].coef.n_elem;
    int nmodels = fitList.size(); 

    mat coefMat(nmodels,ncoef);
    mat seMat(nmodels,ncoef);
    NumericVector rdf(nmodels);
    NumericVector sigSq(nmodels);
    vector<string> ID;
    ID.reserve(nmodels);

    // check optional properties of model
    bool useVCOV = (fitList[0].vcov.n_rows !=0);
    bool useResid = (fitList[0].residuals.n_elem !=0);
    bool useHat = (fitList[0].hatvalues.n_elem !=0);

    mat vcovStacked;
    mat residuals, hatvalues;

    if( useVCOV ){
      // column i stores the vcov matrix for response i
      vcovStacked = mat(pow(ncoef,2), nmodels);
    }

    if( useResid ){
      residuals = mat(fitList[0].residuals.n_elem, nmodels);
    }

    if( useHat ){
      hatvalues = mat(fitList[0].hatvalues.n_elem, nmodels);
    }

    // for each model
    for(int i=0; i< fitList.size(); i++){
      coefMat.row(i) = fitList[i].coef.t();
      seMat.row(i) = fitList[i].se.t();
      sigSq(i) = fitList[i].sigSq;
      rdf(i) = fitList[i].rdf;
      ID.push_back(fitList[i].ID);
      if( useVCOV ) vcovStacked.col(i) = fitList[i].vcov.as_col();
      if( useResid ) residuals.col(i) = fitList[i].residuals;
      if( useHat ) hatvalues.col(i) = fitList[i].hatvalues;
    }

    // Create standard return values
    CharacterVector IDcv(ID.begin(), ID.end());

    sigSq.attr("names") = IDcv;
    rdf.attr("names") = IDcv;

    List lst = List::create(
        Named("coef") = coefMat,
        Named("se") = seMat,
        Named("sigSq") = sigSq,
        Named("rdf") = rdf
      );

    // Insert other return values
    if( useVCOV ) lst["vcovStacked"] = vcovStacked;
    if( useResid ) lst["residuals"] = residuals;
    if( useHat ) lst["hatvalues"] = hatvalues;

    // Set names
    rownames(lst["coef"]) = IDcv;
    rownames(lst["se"]) = IDcv;;
    if( useVCOV )   colnames(lst["vcovStacked"]) = IDcv;
    if( useResid )  colnames(lst["residuals"]) = IDcv;
    if( useHat )    colnames(lst["hatvalues"]) = IDcv;

    return lst;
}






//' Fit series of linear regression models with the same response 
//'
//' Fit regression model \code{y ~ X_design + X_features[,j]} for each feature j
//'
//' @usage lmFitFeatures(y, X_design, X_features, ids, weights, detail = 0L, preprojection = TRUE, nthreads = 1L)
//'
//' @param y response vector
//' @param X_design design matrix shared across all model
//' @param X_features feature matrix with model j using feature j
//' @param ids array of features ids storing identifier for each column in X_features
//' @param weights sample-level weights
//' @param detail level of model detail returned, with LOW = 0, MEDIUM = 1, HIGH = 2. LOW (beta, se, sigSq, rdf), MEDIUM (vcov), HIGH (residuals), MOST (hatvalues)
//' @param preprojection default TRUE. Use preproject of design matrix to accelerate calculations
//' @param nthreads number of threads.  Each model is fit in serial, analysis is parallelized across features
//' 
//' @return List of parameter estimates with entries \code{coef},  \code{se}, \code{sigSq}, \code{rdf} and other depending on \code{detail}
//' @name lmFitFeatures
//' 
//' @examples
//' n = 100 # number of samples
//' p = 10  # number of features
//' nc = 3  # number shared covariates
//' set.seed(1)
//' y = rnorm(n)
//' X = matrix(rnorm(n*p), n, p)
//' colnames(X) = seq(p)
//' X_design = matrix(rnorm(n*nc), n,nc)
//' w = seq(n)
//' w = w / mean(w)
//' 
//' # fit regressions with model j including X[,j]
//' fit = lmFitFeatures(y, X_design, X, colnames(X), w)
//' 
//' # examine results
//' lapply(fit, head, 2)
//' 
//' @export
//
// [[Rcpp::export("lmFitFeatures")]]
List lmFitFeatures_export(const arma::vec &y, const arma::mat &X_design, const arma::mat &X_features, const vector<string> &ids, const arma::vec &weights, const int detail = 0, const bool &preprojection = true, const int &nthreads = 1){

  if( detail > 3) stop("detail > 3 not defined");

  ModelDetail md = static_cast<ModelDetail>(detail);

  if( preprojection && md == MOST ){
    stop("Cannot compute hatvalues with pre-projection");
  }

  vector<ModelFit> fitList = lmFitFeatures(y, X_design, X_features, ids, weights, md, preprojection, nthreads);

  return toList(fitList);
}


//' Fit series of linear regression models to multiple responses with shared design matrix  
//'
//' Fit regression model \code{Y[,j] ~ X_design} for each feature j
//'
//' @param Y matrix of responses as columns
//' @param X design matrix
//' @param ids vector<string> storing identifier for each column in Y
//' @param Weights matrix sample-level weights the same dimension as Y
//' @param detail level of model detail returned, with LOW = 0, MEDIUM = 1, HIGH = 2. LOW (\code{beta}, \code{se}, \code{sigSq}, \code{rdf}), MEDIUM (\code{vcov}), HIGH (\code{residuals}), MOST (\code{hatvalues})
//' @param nthreads number of threads.  Each model is fit in serial, analysis is parallelized across responses.
//'  
//' @details Since the weights vary for each response, each model is computed separately without recycling precomputed values
//' 
//' @return List of parameter estimates with entries \code{coef},  \code{se}, \code{sigSq}, \code{rdf} and other depending on \code{detail}
//' 
//' @name lmFitResponses
//' @examples
//' n = 100
//' m = 5
//' nc = 2
//' set.seed(1)
//' Y = matrix(rnorm(n*m), n, m)
//' X = matrix(rnorm(n*nc), n,nc)
//' colnames(Y) = seq(m)
//' W = matrix(runif(n*m), n,m) 
//' 
//' # fit regressions with model j using Y[,j] as a response
//' fit = lmFitResponses(Y, X, colnames(Y), W)
//' 
//' # examine results
//' lapply(fit, head, 2)
//' 
//' @export
// [[Rcpp::export("lmFitResponses")]]
List lmFitResponses_export(const arma::mat &Y, const arma::mat &X, const vector<string> &ids, const arma::mat &Weights, const int detail = 0,const int &nthreads = 1){

  if( detail > 3) stop("detail > 3 not defined");

  ModelDetail md = static_cast<ModelDetail>(detail);

  vector<ModelFit> fitList = lmFitResponses(Y, X, ids, Weights, md, nthreads);

  return toList(fitList);
}

