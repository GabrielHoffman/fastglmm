#include <armadillo>

#include "local_min.h"

#ifndef NB_THETA_ML_H_
#define NB_THETA_ML_H_

using namespace arma;

/** Store number of observations for each count
 */ 
typedef unordered_map<long, double> CountTable;

/** Compute number of observations for each count
 */ 
static CountTable CreateLUT( const vec &y, const vec &weights = {}, const int &maxUniqueEntries = 1000){

	vec w(weights);
	if( w.is_empty() ){
		w = vec(y.n_elem, fill::ones);
	}

	// lookup table for suming weights for each y value
	unordered_map<long, double> lut;
	lut.reserve(maxUniqueEntries);

	// for each entry in y,
	// increament count for that index
	for(int i=0; i<y.n_elem; i++){
		// increment count for y[i]
		// ++lut[(long) y[i]];
		
		// increment weights for y[i]
		lut[(long) y[i]] += w[i];

		// if number of entries exceeds maxUniqueEntries
		// clear the table and break
		// to return empty table
		if( lut.size() > maxUniqueEntries){
			lut.clear();
			break;
		}
	}

	return lut;
}


struct nbData {
	vec y;
	vec mu;
	double n;
	vec weights;
	mat X;
	bool doCoxReid;
	CountTable ct;	

	// Constructor
	nbData(){}

	// Constructor
	nbData(const vec &y, const vec &mu, const double &n, const vec &weights, const mat &X, const bool &doCoxReid, const CountTable &ct = {}	):
		y(y), mu(mu), n(n), weights(weights), X(X), doCoxReid(doCoxReid), ct(ct) {}
};

// log-likelihood for NB GLM
static  double nb_ll( 	const vec &y,
				const vec &mu,
				const double &n, 
				const vec &weights, 
				const mat &X,
				const bool &doCoxReid,
				const CountTable &ct,	
				const double &theta){

	vec y_theta = y + theta;

	double res;

	if( ct.size() == 0){	
		// vanilla
		vec value = lgamma(y_theta) - lgamma(theta) + theta*log(theta) - y_theta%log(mu + theta);
		res = dot(value, weights) / n;
	}else{
		// count table
		double sum_lgamma_y_theta = 0;
		for(auto &it: ct){
			sum_lgamma_y_theta += it.second * lgamma(it.first + theta);
		}

		double sum_w = sum(weights);

		res = sum_lgamma_y_theta / n - 
				sum_w*lgamma(theta) / n + 
				sum_w*theta*log(theta)/n  - 
				dot(weights,y_theta%log(mu + theta))/n;
	}

	// Doesn't change with theta,
	// so don't need to evaluate it
	// value += - lgamma(y+1) + y%log(mu);

	double cr = 0;
	if( doCoxReid ){
		// Compute logDet of information matrix
		// use sample weights here
		// use 1/theta here since notation is different
		// than in glmGammaPoi
		vec w = weights / (1.0/mu + 1.0/theta);
		double ld;
		bool success = log_det_sympd(ld, X.t() * (X.each_col()%w));
		cr = -0.5 * ld * 0.99; 
		if( ! success ) cr = 0;
	}

	res += cr / n;

	return res;
}

// define functions to be minimized
static inline double nb_ll( double theta_log, void *arg){

	// cast void pointer to nbData
  auto *data = (nbData *) arg;

  // value function to be minimized
  return -1.0*nb_ll(data->y, data->mu, data->n, data->weights, data->X, data->doCoxReid, data->ct, exp(theta_log));
}




/*
MASS::theta.ml() and glmGamPois() maximize the likelihood function using Fishers method with the score and info functions.  That has the advantage of 1) using a starting value from the last call, 2) computes a look up table to compute lgamma(y + theta) since y has a small number of unique values, and 3) gives an estimate of the standard error of theta
*/

/** Estimate dispersion parameter for negative binomial distribution.  Based on \code{MASS::theta.ml()}, but directly maximizes log-likelihood of negative binomial with respect to theta using a line search.
 */ 
static double nb_theta_ml(	
				const vec &y,
				const vec &mu,
				const double &n, 
				const vec &weights = {},		
				const mat &X = {},	
				const bool doCoxReid = true,
				const CountTable &ct = {},	
				const double &left = -5,				
				const double &right = 20,				
				const double &tol = 0.0001220703){

	double theta_log;
	int iter = 0;

	// Setup code to minimize function

	// initialize function
	funcStruct F; 
	F.function = nb_ll;

	vec w(weights);

	if( w.is_empty() ){
		w = vec(y.n_elem, fill::ones);
	}
	nbData d(y, mu, n, w, X, doCoxReid, ct); 
	F.params = &d;

	// maximize log-likelihood
	// theta_log is returned by reference
	local_min(left, right, tol, &F, theta_log, iter);

	return exp(theta_log);
}	

#endif
