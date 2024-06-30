#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifndef THETA_ML_H_
#define THETA_ML_H_

using namespace Rcpp;

// June 29, 2024
// MASS::theta.ml() is 2x faster
// only reason to use this is if I need it in C++ 
// under the hood

// Rcpp armadillo implementation of MASS::theta.ml()
// https://github.com/cran/MASS/blob/38faad6b25fa2b2419b1cd06c78942fa7f5fb7a1/R/negbin.R#L328

// score function 
double score( 	const int &n, 
				const double &th, 
				const NumericVector &mu,
				const NumericVector &y,
				const NumericVector &w){
	double value = sum(w*(digamma(th + y) - digamma(NumericVector(y.length(),th)) + log(th) + 1.0 - log(th + mu) - (y + th)/(mu + th)));

	return value;
}

// information function 
double info( 	const int &n, 
				const double &th, 
				const NumericVector &mu,
				const NumericVector &y,
				const NumericVector &w){

	double value = sum(w*( - Rcpp::trigamma(th + y) + Rcpp::trigamma(NumericVector(y.length(),th)) - 1.0/th + 2.0/(mu + th) - (y + th)/pow(mu + th, 2)));
	return value;
 }

// [[Rcpp::export("theta_ml")]]
List theta_ml(	const NumericVector &y,
				const NumericVector &mu,
				const int &n, 
				const NumericVector &weights,
				const int &limit,
				const double &eps){

	double t0 = ((double) n) / sum(weights*pow(y/mu - 1,2));
	int it = 0;
    double del = 1;
    double i;

    while( (it < limit-1) & (abs(del) > eps) ){
    	t0 = abs(t0);
    	i = info(n, t0, mu, y, weights);
        del = score(n, t0, mu, y, weights) / i;
        t0 += del;
        it++;
    }

	return List::create( 
	        Named("theta") = t0, 
	        Named("se") = sqrt(1/i));
}


#endif