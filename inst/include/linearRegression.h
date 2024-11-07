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

// if -D ARMA, use plain armadillo library
#ifdef ARMA
#include <armadillo>
#else
#include <RcppArmadillo.h>
#endif

// DISABLE warning: solve(): system is singular
#define ARMA_WARN_LEVEL 1

using namespace arma;
using namespace std;

namespace fastlmmLib {

/** Store results from fitting linear regression model
 */ 
struct ModelFit {
	vec coef;
	vec se;
	double df;
	string ID;
	vec hatvalues;
	ModelFit() {}
	ModelFit( const vec &coef, const vec &se, const double &df) : 
		coef(coef), se(se), df(df) {}
	ModelFit( const vec &coef, const vec &se, const double &df, const vec &hatvalues) : 
		coef(coef), se(se), df(df), hatvalues(hatvalues) {}
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
 * @param hat bool indicating if hatvalues should be returned
 * @param rdf_offset degrees of freedom to remove due to pre-projection
 * @param work LMWork workspace to store intermediate results
 * 
/ adapted from https://github.com/RcppCore/RcppArmadillo/blob/master/src/fastLm.cpp
/ https://genomicsclass.github.io/book/pages/qr_and_regression.html
*/
static ModelFit lm(const arma::mat& X, const arma::colvec& y, const bool &hat = false, const double &rdf_offset = 0, LMWork *work = nullptr) {

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
   
	ModelFit fit;
	if( hat ){
		vec hatvalues = diagvec(work->Q * trans(work->Q));
		fit = ModelFit( beta, stderr, rdf, hatvalues);
	}else{
		fit = ModelFit( beta, stderr, rdf);
	}

	// free work if allocated in this function
    if( alloc_local) delete work;

	return fit;
}

static tuple<vec, mat> preprojection(const vec &y, const mat &X_design, const mat &X_features, const vec &weights = {}){

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
	mat X_features_wsqrt = X_features.each_col() % wsqrt;

	// Use QR decomp of X_design, and recycle pre-computed values
	mat Q, R;
	qr_econ(Q, R, X_design_wsqrt);

    // back solve
	vec beta = solve(R, trans(Q) * y_wsqrt);
	vec y_proj = y_wsqrt - X_design_wsqrt * beta;

    // back solve
	mat gamma = solve(R, trans(Q) * X_features_wsqrt);
	mat X_proj = X_features_wsqrt - X_design_wsqrt * gamma;

	return {y_proj, X_proj};
}



/** Weighted linear regression model using QR decomposition.  Recycles QR and memory
 * 
 * @param X design matrix
 * @param y response vector
 * @param weights vector of weights
 * @param hat bool indicating if hatvalues should be returned
 * @param rdf_offset degrees of freedom to remove due to pre-projection
 * @param work LMWork workspace to store intermediate results
 * 
 Scale y and X by sqrt(w / mean(w)) and then call lm()
*/
static ModelFit wlm(const arma::mat& X, const arma::colvec& y, const arma::colvec& w = {}, const bool &hat = false, const double &rdf_offset = 0, LMWork *work = nullptr) {

	ModelFit fit;

	if( w.is_empty() ){
		fit = lm( X, y, hat, rdf_offset, work );
	}else{
		arma::colvec wsqrt = sqrt(w / mean(w));
		fit = lm( X.each_col() % wsqrt, y % wsqrt, hat, rdf_offset, work );
	}

	return fit;
}

// include weights here
// include preprojection

/** Fit series of linear regression model  
 * 
 * @param y response vector
 * @param X_design design matrix
 * @param X_features design matrix with additional features to be fit one at a time
 * @param ids vector<string> storing identifier for each column in X_features
 * @param nthreads number of threads.  Each model is fit in serial, analysis is parallelize across features
 * 
*/
vector<ModelFit> lmFitFeatures(const arma::vec &y, const arma::mat &X_design, const arma::mat &X_features, const vector<string> &ids, const arma::vec &weights = {}, const int &nthreads = 1){

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
			ModelFit fit = wlm(X, y, weights, false, 0, work);
			fit.ID = ids[j];

			// save result to list
			fitList.at(j) =  fit;
		}  
		delete work;
	}  

	return fitList;
}



/** Fit series of linear regression models using pre-projection  
 * 
 * @param y response vector
 * @param X_design design matrix
 * @param X_features design matrix with additional features to be fit one at a time
 * @param ids vector<string> storing identifier for each column in X_features
 * @param nthreads number of threads.  Each model is fit in serial, analysis is parallelize across features
 * 
*/
vector<ModelFit> lmFitFeatures_preproj(const arma::vec &y, const arma::mat &X_design, const arma::mat &X_features, const vector<string> &ids, const arma::vec &weights = {}, const int &nthreads = 1){

	int n_covs = X_design.n_cols;

	vector<ModelFit> fitList(X_features.n_cols, ModelFit());
	bool hat = false;  // don't compute hatvalues

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

	#pragma omp parallel
	{
		LMWork *work = new LMWork();

		// iterate through responses 
		#pragma omp for		 
		for(int j=0; j<X_proj.n_cols; j++){
			// linear regression		
			ModelFit fit = lm(X_proj.col(j), y_proj, hat, rdf_offset, work );
			fit.ID = ids[j];

			// save result to list
			fitList.at(j) =  fit;
		}  
		delete work;
	}  

	return fitList;
}




// include weights here

vector<ModelFit> lmFitResponse(const arma::mat &Y, const arma::mat &X, const vector<string> &ids, const arma::mat &Weights, const int &nthreads = 1){

    vector<ModelFit> fitList(Y.n_cols, ModelFit());

    #ifdef _OPENMP 
        // set threads
        omp_set_num_threads(nthreads);
        // disable nested parallelism
        omp_set_max_active_levels(1);
    #endif

    #pragma omp parallel
    {
        // iterate through responses 
        #pragma omp for      
        for(int j=0; j<Y.n_cols; j++){

            // linear regression        
            ModelFit fit = wlm(X, Y.col(j), Weights.col(j));
			fit.ID = ids[j];

            // save result to list
            fitList.at(j) =  fit;
        }  
    }  

    return fitList;
}
}



#endif