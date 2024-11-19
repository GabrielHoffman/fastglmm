/***********************************************************************
 * @file		linearRegression.h
 * @author		Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		Evaluate linear regression with Armadillo library
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/


#ifndef LINEAR_REGRESSION_H_
#define LINEAR_REGRESSION_H_

#include <tuple>
#include <type_traits>
#include "misc.h"

// DISABLE warning: solve(): system is singular
#define ARMA_WARN_LEVEL 1
// Important: this definition ensures Armadillo enables SuperLU
// #define ARMA_USE_SUPERLU 1

// if -D ARMA, use plain armadillo library
#ifdef ARMA
#include <armadillo>
#else
#include <RcppArmadillo.h>
#endif

using namespace arma;
using namespace std;

namespace fastlmmLib {

// Specify level of model detail to return from regression fit
typedef enum {
    LOW, 	// baseline parameters: beta, se, sigSq, rdf
    MEDIUM, // vcov
    HIGH,   // residuals
    MOST    // hatvalues
} ModelDetail;

/** Store results from fitting linear regression model
 */ 
struct ModelFit {
	vec coef;
	vec se;
	double sigSq;
	double rdf;
	string ID;
	mat vcov;
	vec residuals;
	vec hatvalues;

	ModelFit() {}

	// LOW
	ModelFit( const vec &coef, const vec &se, const double & sigSq, const double &rdf) : 
		coef(coef), se(se), sigSq(sigSq), rdf(rdf) {}

	// MEDIUM
	ModelFit( const vec &coef, const vec &se, const double & sigSq, const double &rdf, const mat & vcov) :
		coef(coef), se(se), sigSq(sigSq), rdf(rdf), vcov(vcov) {}

	// HIGH
	ModelFit( const vec &coef, const vec &se, const double & sigSq, const double &rdf, const mat & vcov, const vec &residuals) : 
		coef(coef), se(se), sigSq(sigSq), rdf(rdf), vcov(vcov), residuals(residuals) {}

	// MOST
	ModelFit( const vec &coef, const vec &se, const double & sigSq, const double &rdf, const mat & vcov, const vec &residuals, const vec &hatvalues) : 
		coef(coef), se(se), sigSq(sigSq), rdf(rdf), vcov(vcov), residuals(residuals), hatvalues(hatvalues) {}
};

/** Workspace for lm() and wlm()
 */ 
struct LMWork {
    mat Q, R, V;
    vec residuals;
    LMWork() {}
};


/** Linear regression model using QR decomposition.  Recycles QR and memory
 * 
 * @param X design matrix
 * @param y response vector
 * @param md return model with specified level of detail. LOW (beta, se, sigSq, rdf), MEDIUM (vcov), HIGH (residuals), MOST (hatvalues)
 * @param rdf_offset degrees of freedom to remove due to pre-projection
 * @param work LMWork workspace to store intermediate results
 * 
/ adapted from https://github.com/RcppCore/RcppArmadillo/blob/master/src/fastLm.cpp
/ https://genomicsclass.github.io/book/pages/qr_and_regression.html
*/
static ModelFit lm(const arma::mat& X, const arma::colvec& y, const ModelDetail md = LOW, const double &rdf_offset = 0, LMWork *work = nullptr) {

	int n = X.n_rows, k = X.n_cols;

	// allocate work, if not already alloc'd
    bool alloc_local = false;
    if( work == nullptr ){
    	alloc_local = true;
        work = new LMWork();
    }

    // QR decomp
	qr_econ(work->Q, work->R, X);

    // back solve
	vec beta = solve(work->R, trans(work->Q) * y);

    // residuals
	work->residuals = y - X*beta;

    // std.errors of coefficients
    double rdf = n - k - rdf_offset;
	double s2 = dot(work->residuals, work->residuals) / rdf; 

    // V = solve(t(R)*R)
    bool success = inv_sympd(work->V, trans(work->R)*work->R);

    vec stderr;
    if( success ){
    	stderr = sqrt(s2 * diagvec(work->V));
    }else{
        stderr = vec(k, fill::value(datum::nan));
        beta.fill(datum::nan);
    }
   
   	// return results with specified level of detail
	ModelFit fit;
    switch( md ){
	    case LOW:
		fit = ModelFit( beta, stderr, s2, rdf);
		break;

		case MEDIUM:
		fit = ModelFit( beta, stderr, s2, rdf, work->V * s2);
		break;

		case HIGH:
		fit = ModelFit( beta, stderr, s2, rdf, work->V * s2, work->residuals);
		break;

		case MOST:
		vec hatvalues = diagvec(work->Q * trans(work->Q));
		fit = ModelFit( beta, stderr, s2, rdf, work->V * s2, work->residuals, hatvalues);
		break;
	}

	// free work if allocated in this function
    if( alloc_local) delete work;

	return fit;
}

/** Regress out covariates using pre-projection
 * 
 * @param y response vector
 * @param X_design design matrix, dense
 * @param X_features matrix with additional features, mat or sp_mat
 * @param weights sample-level weights
 * 
*/
template <typename T>
static tuple<vec, T> preprojection(const vec &y, const mat &X_design, const T &X_features, const vec &weights = {}){

	// naive calculation
	// vec y_proj = y - X_design * inv(trans(X_design) * X_design) * trans(X_design) * y;
	// mat X_proj = X_features - X_design * inv(trans(X_design) * X_design) * trans(X_design) * X_features;

	// if weights is empty, set w to ones
	vec w;
	if( ! weights.is_empty() ){
		w = weights;
	}else{
		w = vec(y.n_elem).ones();
	}
	arma::colvec wsqrt = sqrt(w / mean(w));

	// apply weights in computation to y, X_design, and X_features
	mat X_design_wsqrt = X_design.each_col() % wsqrt;
	vec y_wsqrt = y % wsqrt;
	T X_features_wsqrt = scaleEachCol(X_features, wsqrt);

	// Use QR decomp of X_design, and recycle pre-computed values
	mat Q, R;
	qr_econ(Q, R, X_design_wsqrt);

    // back solve
	vec beta = solve(R, trans(Q) * y_wsqrt);
	vec y_proj = y_wsqrt - X_design_wsqrt * beta;

    // back solve
    // use constructor T() to subtract matricies of the same type
	mat gamma = solve(R, trans(Q) * X_features_wsqrt);
	T X_proj;

	// cast X_design_wsqrt * gamma to type T if needed
	if( is_same_v<decltype(X_features_wsqrt), decltype(X_design_wsqrt)> ){
		X_proj = X_features_wsqrt - X_design_wsqrt * gamma;
	}else{
		X_proj = X_features_wsqrt - T(X_design_wsqrt * gamma);
	}

	return {y_proj, X_proj};
}


/** Regress out covariates using pre-projection
 * 
 * @param y response vector
 * @param X_design design matrix, sparse
 * @param X_features matrix with additional features, mat or sp_mat 
 * @param weights sample-level weights
*/
template <typename T>
static tuple<vec, T> preprojection(const vec &y, const sp_mat &X_design, const T &X_features, const vec &weights = {}){

	// naive calculation
	// vec y_proj = y - X_design * inv(trans(X_design) * X_design) * trans(X_design) * y;
	// mat X_proj = X_features - X_design * inv(trans(X_design) * X_design) * trans(X_design) * X_features;

	// if weights is empty, set w to ones
	vec w;
	if( ! weights.is_empty() ){
		w = weights;
	}else{
		w = vec(y.n_elem).ones();
	}
	arma::colvec wsqrt = sqrt(w / mean(w));

	// apply weights in computation to y, X_design, and X_features
	sp_mat X_design_wsqrt = scaleEachCol(X_design, wsqrt);
	vec y_wsqrt = y % wsqrt;
	T X_features_wsqrt = scaleEachCol(X_features, wsqrt);

	// recycle sparse crossprod
	// spsolve uses lapack after converting V to dense matrix
	// minimal penalty practical V dimensions
	bool success;
	vec beta;
	mat gamma;

	sp_mat V = trans(X_design_wsqrt) * X_design_wsqrt;
	success = spsolve(beta, V, vec(trans(X_design_wsqrt) * y_wsqrt), "lapack");
	if( ! success ) beta.fill(datum::nan);
	vec y_proj = y_wsqrt - X_design_wsqrt * beta;

	success = spsolve(gamma, V, mat(trans(X_design_wsqrt) * X_features_wsqrt), "lapack");
	if( ! success ) gamma.fill(datum::nan);

	T X_proj;
	// cast X_design_wsqrt * gamma to type T if needed
	if( is_same_v<decltype(X_features_wsqrt), decltype(X_design_wsqrt)> ){
		X_proj = X_features_wsqrt - X_design_wsqrt * gamma;
	}else{
		X_proj = X_features_wsqrt - T(X_design_wsqrt * gamma);
	}

	// Use LU decomp of X_design
	// requires SuperLU library
	//  -L/opt/homebrew/Cellar/superlu/7.0.0/lib -lsuperlu
	// but not portable
	// spsolve_factoriser SF;

	// bool status = SF.factorise(X_design_wsqrt);

	// vec beta;
	// SF.solve(beta, y_wsqrt);
	// vec y_proj = y_wsqrt - X_design_wsqrt * beta;

	// vec gamma;
	// SF.solve(gamma, X_features_wsqrt);
	// T X_proj = X_features_wsqrt - T(X_design_wsqrt * gamma);

	return {y_proj, X_proj};
}



/** Weighted linear regression model using QR decomposition.  Recycles QR and memory
 * 
 * @param X design matrix
 * @param y response vector
 * @param weights vector of weights
 * @param md return model with specified level of detail. LOW (beta, se, sigSq, rdf), MEDIUM (vcov), HIGH (residuals), MOST (hatvalues)
 * @param rdf_offset degrees of freedom to remove due to pre-projection
 * @param work LMWork workspace to store intermediate results
 * 
 Scale y and X by sqrt(w / mean(w)) and then call lm()
*/
static ModelFit wlm(const arma::mat& X, const arma::colvec& y, const arma::colvec& w = {}, const ModelDetail md = LOW, const double &rdf_offset = 0, LMWork *work = nullptr) {

	ModelFit fit;

	if( w.is_empty() ){
		fit = lm( X, y, md, rdf_offset, work );
	}else{
		arma::colvec wsqrt = sqrt(w / mean(w));
		fit = lm( X.each_col() % wsqrt, y % wsqrt, md, rdf_offset, work );

		if( md >= HIGH){
            // Rescale residuals by weights afterward
            //  since input X and y are scaled before lm()
            fit.residuals /= wsqrt;
        }
	}

	return fit;
}

/** Fit series of linear regression models  
 * 
 * @param y response vector
 * @param X_design design matrix
 * @param X_features matrix with additional features to be fit one at a time
 * @param ids vector<string> storing identifier for each column in X_features
 * @param weights sample-level weights
 * @param md return model with specified level of detail. LOW (beta, se, sigSq, rdf), MEDIUM (vcov), HIGH (residuals), MOST (hatvalues)
 * @param nthreads number of threads.  Each model is fit in serial, analysis is parallelized across features
 * 
*/
vector<ModelFit> lmFitFeatures_standard(const arma::vec &y, const arma::mat &X_design, const arma::mat &X_features, const vector<string> &ids, const arma::vec &weights = {}, const ModelDetail md = LOW, const int &nthreads = 1){

	int n_covs = X_design.n_cols;

	vector<ModelFit> fitList(X_features.n_cols, ModelFit());
	
	#ifdef _OPENMP 
		// set threads
		omp_set_num_threads(nthreads);
		// disable nested parallelism
		omp_set_max_active_levels(1);
	#endif

	#pragma omp parallel
	{
		// create design matrix with jth feature in the last column
		// X = cbind(X_design, X_features[,0])
		arma::mat X(X_design);
		X.insert_cols(n_covs, X_features.col(0));

		LMWork *work = new LMWork();

		// iterate through responses 
		#pragma omp for		 
		for(int j=0; j<X_features.n_cols; j++){
			// Create design matrix with intercept as first column
			X.col(n_covs) = X_features.col(j);

			// linear regression		
			ModelFit fit = wlm(X, y, weights, md, 0, work);
			fit.ID = ids[j];

			// save result to list
			fitList.at(j) =  fit;
		}  
		delete work;
	}  

	return fitList;
}



/** Fit series of linear regression models using pre-projection.  
 * 
 * @param y response vector
 * @param X_design design matrix, mat or sp_mat
 * @param X_features design matrix, mat or sp_mat, with additional features to be fit one at a time
 * @param ids vector<string> storing identifier for each column in X_features
 * @param weights sample-level weights
 * @param md return model with specified level of detail. LOW (beta, se, sigSq, rdf), MEDIUM (vcov), HIGH (residuals), MOST (hatvalues)
 * @param nthreads number of threads.  Each model is fit in serial, analysis is parallelized across features
 * 
*/
template <typename T1, typename T2>
vector<ModelFit> lmFitFeatures_preproj(const arma::vec &y, const T1 &X_design, const T2 &X_features, const vector<string> &ids, const arma::vec &weights = {}, const ModelDetail md = LOW, const int &nthreads = 1){

	vector<ModelFit> fitList(X_features.n_cols, ModelFit());
	
	#ifdef _OPENMP 
		// set threads
		omp_set_num_threads(nthreads);
		// disable nested parallelism
		omp_set_max_active_levels(1);
	#endif

	// pre-projection to regress out the covariates first
	// set rdf_offset to the number of covariates projected out
	auto [y_proj, X_proj] = preprojection(y, X_design, X_features, weights);
	double rdf_offset = X_design.n_cols;

	vec wsqrt = sqrt(weights);

	#pragma omp parallel
	{
		LMWork *work = new LMWork();

		// iterate through responses 
		#pragma omp for		 
		for(int j=0; j<X_proj.n_cols; j++){
			// linear regression		
			ModelFit fit = lm(X_proj.col(j), y_proj, md, rdf_offset, work );

			fit.ID = ids[j];

			if( md >= HIGH){
	            // Rescale residuals by weights afterward
	            //  since input X and y are scaled before lm()
	            fit.residuals /= wsqrt;
	        }

			// save result to list
			fitList.at(j) =  fit;
		}  
		delete work;
	}  

	return fitList;
}

/** Fit series of linear regression models using pre-projection.  
 * 
 * @param y response vector
 * @param X_design design matrix, mat or sp_mat
 * @param X_features design matrix, mat or sp_mat, with additional features to be fit one at a time
 * @param ids vector<string> storing identifier for each column in X_features
 * @param weights sample-level weights
 * @param md return model with specified level of detail. LOW (beta, se, sigSq, rdf), MEDIUM (vcov), HIGH (residuals), MOST (hatvalues)
 * @param preprojection default true. Use preproject of design matrix to accelerate calculations
 * @param nthreads number of threads.  Each model is fit in serial, analysis is parallelized across features
 * 
*/
template <typename T1, typename T2>
vector<ModelFit> lmFitFeatures(const arma::vec &y, const T1 &X_design, const T2 &X_features, const vector<string> &ids, const arma::vec &weights = {}, const ModelDetail md = LOW, const bool &preprojection = true, const int &nthreads = 1){

	vector<ModelFit> fitList;

	if( preprojection ){
		// supports mat and sp_mat
		fitList = lmFitFeatures_preproj(y, X_design, X_features, ids, weights, md, nthreads);
	}else{
		// only supports mat 
		fitList = lmFitFeatures_standard(y, mat(X_design), mat(X_features), ids, weights, md, nthreads);
	}

	return fitList;
}






/** Fit series of linear regression models to multiple responses with shared design matrix  
 * 
 * @param Y matrix of responses as columns
 * @param X design matrix
 * @param ids vector<string> storing identifier for each column in Y
 * @param Weights matrix sample-level weights the same dimension as Y
 * @param md return model with specified level of detail. LOW (beta, se, sigSq, rdf), MEDIUM (vcov), HIGH (residuals), MOST (hatvalues)
 * @param nthreads number of threads.  Each model is fit in serial, analysis is parallelized across responses.
 * 
 * Since the weights vary for each response, each model is computed separately without recycling precomputed values
*/
vector<ModelFit> lmFitResponses(const arma::mat &Y, const arma::mat &X, const vector<string> &ids, const arma::mat &Weights, const ModelDetail md = LOW, const int &nthreads = 1){

    vector<ModelFit> fitList(Y.n_cols, ModelFit());

    #ifdef _OPENMP 
        // set threads
        omp_set_num_threads(nthreads);
        // disable nested parallelism
        omp_set_max_active_levels(1);
    #endif

    arma::mat Wsqrt = sqrt(Weights);
    arma::mat Yw = Y % Wsqrt;

    #pragma omp parallel
    {
        // iterate through responses 
        #pragma omp for      
        for(int j=0; j<Y.n_cols; j++){

            // linear regression        
            // ModelFit fit = wlm(X, Y.col(j), Weights.col(j));
            ModelFit fit = lm(X.each_col() % Wsqrt.col(j), Yw.col(j), md);

			fit.ID = ids[j];

			if( md >= HIGH){
	            // Rescale residuals by weights afterward
	            //  since input X and y are scaled before lm()
	            fit.residuals /= Wsqrt.col(j);
        	}

            // save result to list
            fitList.at(j) =  fit;
        }  
    }  

    return fitList;
}

}



#endif