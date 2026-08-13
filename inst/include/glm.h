/***************************************************************
 * @file		glm.h
 * @author	Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		Fit generalizd linear models
 * Copyright (C) 2024 Gabriel Hoffman	
 **************************************************************/


#ifndef _GLM_H_
#define _GLM_H_

// if -D USE_R, use RcppArmadillo library
#ifdef USE_R
// [[Rcpp::depends(RcppParallel)]]  
// [[Rcpp::depends(RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#else
#include <armadillo>
#include <tbb/tbb.h>
#endif


#include <iostream>

// #include "ModelFit.h"
#include "linearRegression.h"
#include "nb_theta.h"
#include "glm_family.h"

using namespace std;
using namespace arma;

namespace fastglmmLib {

/** Check that response is valid for given family
 */ 
static void checkResponse(const vec &y, const string &family){
	
	shared_ptr<GLMFamily> fam = getGLMFamily( family );

	// Keep unique sorted values that are not NAN
	vec res = unique(omit_nan(y));

	string famStr = fam->family();

	if( famStr == "BinomialLogit" || famStr == "BinomialProbit" ){
		// valid for binary logistic regression
		bool valid_binary = (res.n_elem == 2 && (res[0] == 0 && res[1] == 1));

		// valid for binomial with logit/probit link
		bool valid_beta = (res[0] >= 0 && res[res.n_elem-1] <= 1);

		if( ! (valid_binary || valid_beta) ){
			throw logic_error( "Invalid response for binomial, must be only 0/1 or in interval [0,1]" );
		}
	}
	else if( famStr == "PoissonLog" || famStr == "QuasipoissonLog" || famStr == "NB"){
		// min value must be non-negative
		if( res[0] < 0){
			throw logic_error( "Invalid response for poisson/nb: must be non-negative" );
		}

		// check that values are integers
		if( ! all_integer_valued(res) ){			
			throw logic_error( "Invalid response for poisson/nb: must be integers" );
		}

	}else if( famStr == "QuasibinomialLogit" ){
		if( res[0] < 0 || res[res.n_elem-1] >=1 ){
			throw logic_error( "Invalid response for quasi-binomial, must be between 0 and 1" );
		}
	}
}


/** Workspace for GLM 
 */ 
struct GLMWork : LMWork {	
	vec eta, mu, gprime, z, wsqrt, w;
    GLMWork() {}
};


/** Generalized linear model fit with IRLS 
 * 
 * Adapted from irls_qrnewton() in https://bwlewis.github.io/GLM/ Values are mu are initialized as in R's glm().  work->residuals stores the Pearson residuals.  Calls lm() on weighted working response and design matrix.
 * 
 * @param X design matrix
 * @param y response vector 
 * @param family type of GLM: "gaussian", "gaussian/identity", "binomial/logit", "binomial/probit", "poisson/log", "quasibinomial/logit", "quasipoisson/log", "nb", nb:x" where x is a numeric value of theta
 * @param md level of model detail returned, with LEAST = 0, LOW = 1, MEDIUM = 2, HIGH = 3, MOST = 4, MAX = 5. LEAST (beta), LOW (beta, se, dispersion, rdf), MEDIUM (vcov), HIGH (pearson residuals), MOST (hatvalues), MAX (deviance residuals)
 * @param weights sample-level weights
 * @param offset sample-level offsets
 * @param work GLMWork workspace to store intermediate results
 * @param betaInit initial value for coefficients. If empty, use initialize() methods from GLMFamily
 * @param epsilon stopping criteria for norm between coefficient estimates compare to the previous iteration
 * @param maxit max iterations
 * @param lambda ridge parameter
 */  
static ModelFitGLM GLM(
				const mat& X, 
				const colvec& y, 
				const string &family, 
				const ModelDetail md = LOW, 
				const vec &weights = {}, 
				const vec &offset = {}, 
				GLMWork *work = nullptr, 
				const vec &betaInit = {}, 
				const double &epsilon = 1e-8, 
				const double &maxit = 25,
				const double &lambda = 0){

	shared_ptr<GLMFamily> fam = getGLMFamily( family );

	checkResponse(y, family);

	// allocate work, if not already alloc'd
  bool alloc_local = false;
  if( work == nullptr ){
  	alloc_local = true;
      work = new GLMWork();
  }

	ModelFit fit;
	fit.coef = betaInit;

	// if offset is empty, set to zeros
	vec offset_(offset);
	if( offset.is_empty() ){		
		offset_ = vec(y.n_elem, fill::zeros);
	}

	// if weights is empty, set to ones1
	vec weights_(weights);
	if( weights.is_empty() ){		
		weights_ = vec(y.n_elem, fill::ones);
	}

	int i;
	for(i=0; i<maxit; i++){    	
		if( i == 0 && betaInit.is_empty() ){			
			// initialize mu, essential for Poisson
			// faster convergence than initializing beta to zero
			work->mu 	= fam->initialize(y, weights_);
			work->eta = fam->link(work->mu) ;
		}else{
	    	// linear predictor
			work->eta = X * fit.coef + offset_;
			work->mu 	= fam->linkinv( work->eta );
		}

  	work->gprime= fam->mu_eta( work->eta );
  	work->z  		= (work->eta - offset_) + (y - work->mu) / work->gprime;

  	// if( fam->family() == "NB"){
		// 	// Numerically stable evaluation for NB		
		// 	double theta = fam->getOverdispersion();		
	  //   vec x = work->eta - log(theta);
	  //   work->wsqrt = sqrt(weights_ %(theta * (0.5 * (1.0 + tanh(0.5 * x)))));
  	// }else{
	  	work->wsqrt = work->gprime % sqrt(weights_ / fam->variance( work->mu ));
	  // }

  	// if weights are zero, set them to very small value
		// uvec idx = find(work->wsqrt == 0 && weights_ != 0);
		// work->wsqrt.elem(idx).fill(1e-8);

  	vec beta_prev(fit.coef);

  	// Solve least squares system to get beta
  	fit = lm(scaleEachCol(X, work->wsqrt), work->z % work->wsqrt, LEAST, lambda, 0, work );

  	// if model is singular
  	if( ! fit.success ){
  		break;
  	}

		// stopping criterion
		if( i > 0 && norm(fit.coef - beta_prev) < epsilon ) break;
	}

	// for last beta value
  // linear predictor
	work->eta = X * fit.coef + offset_;
	work->mu = fam->linkinv( work->eta );

  // Solve least squares system, 
  // estimate other parameters based on ModelDetail
  // Estimate dispersion if needed
  // Reduce residual degrees of freedom by the number of 
  // 	entries with zero weights
  double rdf_offset = sum(work->wsqrt == 0);
  fit = lm(scaleEachCol(X, work->wsqrt), work->z % work->wsqrt, md, lambda, rdf_offset, work, fam->estimateDispersion(), true);

  if( md == MAX){
		// compute raw deviance residuals
  	vec dr = fam->dev_resids(y, work->mu, weights_);

  	// transform and store residuals
  	fit.setDevResids( dr, y, work->mu, weights_);
  }

  if( md >= MOST ){
  	fit.setFittedValues( work->mu, weights_ );
  }

  if( md >= HIGH ){
	  // if weight is zero, set residuals to NAN
  	fit.residuals.elem(find(work->wsqrt == 0)).fill(datum::nan);
  }

  if( md >= LOW ){
		// save variance of fitted eta:
		// prediction on latent scale
	  fit.varFitted = wvar( work->eta, weights_);
  }

	// free work if allocated in this function
  if( alloc_local) delete work;

	return ModelFitGLM(fit, family, i);
}






/** Negative binomial fit
 * 
 * Fit Poisson model, then estimate theta parameter of negative binomial.  Coefficients have a small dependence on theta that is apparent with using offset, so need to iteratively estimate coefficients and theta.
 * 
 * @param X design matrix
 * @param y response vector 
 * @param md level of model detail returned, with LEAST = 0, LOW = 1, MEDIUM = 2, HIGH = 3, MOST = 4, MAX = 5. LEAST (beta), LOW (beta, se, dispersion, rdf), MEDIUM (vcov), HIGH (pearson residuals), MOST (hatvalues), MAX (deviance residuals)
 * @param weights sample-level weights
 * @param offset sample-level offsets
 * @param doCoxReid use Cox-Reid adjustment when estimating overdispersion for negative binomial models.
 * @param work GLMWork workspace to store intermediate results
 * @param betaInit initial value for coefficients. If empty, use initialize() methods from GLMFamily
 * @param epsilon stopping criteria for norm between coefficient estimates compare to the previous iteration
 * @param maxit max iterations
  * @param epsilon_nb stopping criteria for norm between coefficient estimates compare to the previous iteration in estimating theta
 * @param maxit_nb max iterations for estimating theta
 * @param lambda ridge parameter
 */  
static ModelFitGLM GLM_NB(	
	const mat& X, 
	const colvec& y, 
	const ModelDetail md = LOW, 
	const vec &weights = {}, 
	const vec &offset = {}, 
	const bool &doCoxReid = true, 
	GLMWork *work = nullptr, 
	const vec &betaInit = {}, 
	const double &epsilon = 1e-8, 
	const double &maxit = 25, 
	const double &epsilon_nb = 1e-4, 
	const double &maxit_nb = 5,
	const double &lambda = 0){

	// allocate work, if not already alloc'd
  bool alloc_local = false;
  if( work == nullptr ){
  	alloc_local = true;
      work = new GLMWork();
  }

	ModelFitGLM fit;
	string family; 
	double theta;

	// Fit Poisson regression to estimate coefficients
	fit = GLM(X, y, "poisson/log", LEAST, weights, offset, work, betaInit, epsilon, maxit, lambda);
	vec beta_prev(fit.coef);

	// Lookup table to speed up estimation of theta
	CountTable ct = CreateLUT(y, weights);

	for(int i=0; i<maxit_nb; i++){

		theta = nb_theta_ml(y, work->mu, y.n_elem, weights, X, doCoxReid, ct);

		// Fit NB regression to estimate coefficients
    beta_prev = fit.coef;
		family = "nb:" + to_string(theta);
		fit = GLM(X, y, family, LEAST, weights, offset, work, fit.coef, epsilon, maxit, lambda);

		// if model is singular
		if( !fit.success ) break;

		// stopping criterion
		if( norm(fit.coef - beta_prev) < epsilon_nb ) break;
	}

	// Estimate parameters based on ModelDetail
	fit = GLM(X, y, family, md, weights, offset, work, fit.coef, epsilon, maxit, lambda);
	fit.theta = theta;

	// Since theta is estimated from the data
	// set the dispersion to be 1
	// So undo scaling by dispersion
	if( md >= MEDIUM ){
		fit.vcov = fit.vcov / fit.dispersion;
		fit.se = sqrt(diagvec(fit.vcov));
	}
	fit.dispersion = 1.0;

	// free work if allocated in this function
  if( alloc_local) delete work;

	return fit;
}







/** Fit series of linear regression models  
 * 
 * @param y response vector
 * @param X_design design matrix
 * @param X_features matrix with additional features to be fit one at a time
 * @param ids vector<string> storing identifier for each column in X_features
 * @param family type of GLM: "gaussian", "gaussian/identity", "binomial/logit", "binomial/probit", "poisson/log", "quasibinomial/logit", "quasipoisson/log", "nb", nb:x" where x is a numeric value of theta
 * @param weights sample-level weights
 * @param offsets sample-level offsets
 * @param md level of model detail returned, with LEAST = 0, LOW = 1, MEDIUM = 2, HIGH = 3, MOST = 4, MAX = 5. LEAST (beta), LOW (beta, se, dispersion, rdf), MEDIUM (vcov), HIGH (pearson residuals), MOST (hatvalues), MAX (deviance residuals)
 * @param doCoxReid use Cox-Reid adjustment when estimating overdispersion for negative binomial models.  Default TRUE for less than 100 samples
 * @param shareTheta estimate theta from design matrix, and share across all features instead of re-estimating for each feature
 * @param fastApprox default false.  if true, use pre-projection on the working response from an initial regression fit on only the X_design.  Under the null for X_features, this is a very good approximation and _much_ faster
 * @param nthreads number of threads.  Each model is fit in serial, analysis is parallelized across features
 * @param epsilon tolerance for GLM IRLS 
 * @param maxit max iterations for GLM IRLS
 * @param epsilon_nb tolerance for negative binomial
 * @param maxit_nb max iterations for negative binomial
 * @param lambda ridge parameter
*/
static ModelFitGLMList glmFitFeatures(	
	const arma::vec &y, 
	const arma::mat &X_design, 
	const arma::mat &X_features, 
	const vector<string> &ids,  
	string family, 
	arma::vec weights = {}, 
	const vec &offset = {}, 
	const ModelDetail md = LOW, 
	const bool &doCoxReid = true, 
	const bool &shareTheta = false, 
	const bool &fastApprox = false,
	const int &nthreads = 1, 
	const double &epsilon = 1e-8, 
	const double &maxit = 25, 
	const double &epsilon_nb = 1e-4,
	const double &maxit_nb = 5,
	const double &lambda = 0){

  // standardize weights
  if( ! weights.is_empty() ){
  	weights = weights / mean(weights);
  }

	int n_covs = X_design.n_cols;

	// Estimate coef using only design matrix
	// use to initialize coefficients for each feature
	ModelFitGLM fitInit;
	GLMWork *work = new GLMWork();
	if( family == "nb" ){
		fitInit = GLM_NB(X_design, y, LEAST, weights, offset, doCoxReid, work, {}, epsilon, maxit, epsilon_nb, maxit_nb, lambda);

		// if shareTheta, use same theta value across all features
		if( shareTheta ){
			family = "nb:" + to_string(fitInit.theta);
		}
	}else{
  	fitInit = GLM(X_design, y, family, LEAST, weights, offset, work, {}, epsilon, maxit, lambda);
  }

  // get working response
  vec workingResponse(work->z);
  vec workingWeights(square(work->wsqrt));
  workingWeights = workingWeights / mean(workingWeights);
  delete work;

	ModelFitGLMList fitList(X_features.n_cols, ModelFitGLM());

  if( fastApprox ){

  	// Pre-projection on working response
  	// when test of X_features is truely under the null,
  	// approximation is very good
		ModelFitList mfl = lmFitFeatures_preproj(workingResponse, X_design, X_features, ids, workingWeights, md, nthreads);

		for(int i=0; i<mfl.size(); i++){
			fitList.at(i) = ModelFitGLM(mfl[i], family, 1);
		}
  }else{
	
  	// Full fit of each model

    // set betaInit to [fitInit.coef,0]
    vec betaInit(fitInit.coef);
    betaInit.resize(betaInit.n_elem + 1);
    betaInit[betaInit.n_elem] = 0;
    
		// Parallel part using Thread Building Blocks
		tbb::task_arena limited_arena(nthreads);
		limited_arena.execute([&] {
		tbb::parallel_for(
			tbb::blocked_range<int>(0, X_features.n_cols, 100), 
			[&](const tbb::blocked_range<int>& r){ 

			disable_parallel_blas();

			// create design matrix with jth feature in the last column
			// X = cbind(X_design, X_features[,0])
			arma::mat X(X_design);
			X.insert_cols(n_covs, X_features.col(0));

			GLMWork *work = new GLMWork();

			// iterate through features 
	  	for (int j = r.begin(); j != r.end(); ++j) {  
				// Create design matrix with intercept as first column
				X.col(n_covs) = X_features.col(j);

				// GLM regression		
				ModelFitGLM fit;
	      	if( family == "nb" ){
	        	fit = GLM_NB(X, y, md, weights, offset, doCoxReid, work, betaInit, epsilon, maxit, epsilon_nb, maxit_nb, lambda);
	      	}else{
	        	fit = GLM(X, y, family, md, weights, offset, work, betaInit, epsilon, maxit, lambda);
	        }

		    // Save feature ID
				fit.ID = ids[j];

				// save result to list
				fitList.at(j) =  fit;
			}  

			delete work;
		}); }); 
	}

	return fitList;
}





/** Fit series of linear regression models to multiple responses with shared design matrix  
 * 
 * @param Y matrix of responses as columns
 * @param X design matrix
 * @param ids vector<string> storing identifier for each column in Y
 * @param family GLM family as a vector of strings, one entry per response 
 * @param weights sample-level weights the same number of samples as Y
 * @param offset vector of sample-level offset values
 * @param md level of model detail returned, with LEAST = 0, LOW = 1, MEDIUM = 2, HIGH = 3, MOST = 4, MAX = 5. LEAST (beta), LOW (beta, se, dispersion, rdf), MEDIUM (vcov), HIGH (pearson residuals), MOST (hatvalues), MAX (deviance residuals)
 * @param nthreads number of threads.  Each model is fit in serial, analysis is parallelized across responses.
 * @param doCoxReid use Cox-Reid adjustment when estimating overdispersion for negative binomial models.  Default TRUE for less than 100 samples
 * @param nthreads number of threads.  Each model is fit in serial, analysis is parallelized across responses.
 * @param epsilon tolerance for GLM IRLS 
 * @param maxit max iterations for GLM IRLS
 * @param epsilon_nb tolerance for negative binomial
 * @param maxit_nb max iterations for negative binomial
 * @param lambda ridge parameter
 * 
 * Since the weights vary for each response, each model is computed separately without recycling precomputed values
*/
static ModelFitGLMList glmFitResponses(
	const arma::mat &Y, 
	const arma::mat &X, 
	const vector<string> &ids, 
	const vector<string> &family, 
	const arma::vec weights = {}, 
	const vec &offset = {}, 
	const ModelDetail md = LOW, 
	const bool &doCoxReid = true, 
	const int &nthreads = 1, 
	const double &epsilon = 1e-8, 
	const double &maxit = 25, 
	const double &epsilon_nb = 1e-4,
	const double & maxit_nb = 5,
	const double &lambda = 0){

	// standardize weights
	vec w_norm = weights;
	if( ! w_norm.is_empty() ){
		w_norm = w_norm / mean(w_norm);
	}

	ModelFitGLMList fitList(Y.n_cols, ModelFitGLM());

	// find rows in X with NAN values
	uvec idx_drop = rows_with_nan(X);  
	mat X_clean(X);
	X_clean.rows(idx_drop).zeros();

	// Parallel part using Thread Building Blocks
	tbb::task_arena limited_arena(nthreads);
	limited_arena.execute([&] {
	tbb::parallel_for(
		tbb::blocked_range<int>(0, Y.n_cols, 10), 
		[&](const tbb::blocked_range<int>& r){ 

		disable_parallel_blas();

		// local workspace 
		GLMWork *work = new GLMWork();
		vec y, w;
		uvec idx;

		// iterate through responses 
		for (int j = r.begin(); j != r.end(); ++j) {  

			// identify samples with NAN entries
			// set values and weights to zero
			y = Y.col(j);	
			w = w_norm;
			idx = unique(join_cols(find_nan(y), idx_drop));
			y.elem(idx).zeros();
			w.elem(idx).zeros();

			// GLM regression   
			ModelFitGLM fit;
			if( family[j] == "nb" ){
				fit = GLM_NB(X_clean, y, md, w, offset, doCoxReid, work, {}, epsilon, maxit, epsilon_nb, maxit_nb, lambda);
			}else{
				fit = GLM(X_clean, y, family[j], md, w, offset, work, {}, epsilon, maxit, lambda);
			}

			// Save feature ID
			fit.ID = ids[j];

			// return mean of mu for jth response
			fit.mu_mean = mean(work->mu);

			// return mean of response
			fit.y_mean = mean(y);

			// save result to list
			fitList.at(j) = fit;
    }  
		delete work;
	}); });

	return fitList;
	}

/** Fit series of linear regression models to multiple responses with shared design matrix  
 * 
 * @param Y matrix of responses as columns
 * @param X design matrix
 * @param ids vector<string> storing identifier for each column in Y
 * @param family GLM family as a string
 * @param weights sample-level weights the same number of samples as Y
 * @param offset vector of sample-level offset values
 * @param md return model with specified level of detail. LOW (beta, se, dispersion, rdf), MEDIUM (vcov), HIGH (residuals), MOST (hatvalues)
 * @param nthreads number of threads.  Each model is fit in serial, analysis is parallelized across responses.
 * @param doCoxReid use Cox-Reid adjustment when estimating overdispersion for negative binomial models.  
 * @param nthreads number of threads.  Each model is fit in serial, analysis is parallelized across responses.
 * @param epsilon tolerance for GLM IRLS 
 * @param maxit max iterations for GLM IRLS
 * @param epsilon_nb tolerance for negative binomial
 * @param maxit_nb max iterations for negative binomial
 * @param lambda ridge parameter
 * 
 * Since the weights vary for each response, each model is computed separately without recycling precomputed values
*/
static ModelFitGLMList glmFitResponses(
	const arma::mat &Y, 
	const arma::mat &X, 
	const vector<string> &ids, 
	const string &family, 
	const arma::vec &weights = {}, 
	const vec &offset = {}, 
	const ModelDetail md = LOW, 
	const bool &doCoxReid = true, 
	const int &nthreads = 1, 
	const double &epsilon = 1e-8, 
	const double &maxit = 25, 
	const double &epsilon_nb = 1e-4, 
	const double & maxit_nb = 5,
	const double &lambda = 0){

	// all responses analyzed with same family value
	vector<string> famVec(Y.n_cols, family);

	return glmFitResponses( Y, X, ids, famVec, weights, offset, md, doCoxReid, nthreads, epsilon, maxit, epsilon_nb, maxit_nb, lambda);
}


};

#endif