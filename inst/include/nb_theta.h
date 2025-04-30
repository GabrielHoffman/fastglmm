#include <armadillo>

#include <local_min.h>

#ifndef NB_THETA_ML_H_
#define NB_THETA_ML_H_

using namespace arma;

struct nbData {
	vec y;
	vec mu;
	double n;
	vec weights;

	// Constructor
	nbData(vec y, vec mu, double n, vec weights):
		y(y), mu(mu), n(n), weights(weights) {}
};

// log-likelihood for NB GLM
double nb_ll( const vec &y,
							const vec &mu,
							const double &n, 
							const vec &weights, 
							const double &theta){

	// compute logLik for each sample
	vec tmp = y + theta;
	vec value = lgamma(tmp) - lgamma(theta) - lgamma(y+1) + y%log(mu) + theta * log(theta) - tmp%log(mu + theta);

	// sum weighted logLik values
	return dot(value, weights) / n;
}


// define functions to be minimized
static inline double nb_ll( double theta_log, void *arg){

	// cast void pointer to nbData
  auto *data = (nbData *) arg;

  // value function to be minimized
  return -1.0*nb_ll(data->y, data->mu, data->n, data->weights, exp(theta_log));
}




double nb_theta_ml(	
				const vec &y,
				const vec &mu,
				const double &n, 
				const vec &weights,				
				const double &left = -5,				
				const double &right = 20,				
				const double &tol = 0.0001220703){

	double theta_log;
  int iter = 0;

  // Setup code to minimize function

	// initialize function
	funcStruct F; 
	F.function = nb_ll;

	nbData d = nbData(y, mu, n, weights);
	F.params = &d;

	// maximize log-likelihood
	// theta_log is returned by reference
	local_min(left, right, tol, &F, theta_log, iter);

	return exp(theta_log);
}	









#endif