/***********************************************************************
 * @file		ModelFit.h
 * @author		Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		Store parameters from model fit
 * Copyright (C) 2024 Gabriel Hoffman
 ***********************************************************************/

#ifndef MODEL_FIT_H_
#define MODEL_FIT_H_

#include <vector>
#include <string>

// from fastlmm
#include <misc.h>

using namespace arma;
using namespace std;

namespace fastlmmLib {

// Specify level of model detail to return from regression fit
typedef enum {
	LEAST,  // just beta
    LOW, 	// baseline parameters: beta, se, dispersion, rdf
    MEDIUM, // vcov
    HIGH,   // pearson residuals
    MOST,   // hatvalues, fitted.values
    MAX     // deviance residuals
} ModelDetail;

/** Store results from fitting linear regression model
 */ 
class ModelFit {
	public: 
	bool success;
	vec coef;
	vec se;
	double dispersion;
	double rdf;
	string ID;
	mat vcov;
	vec residuals;
	vec hatvalues;
	vec mu;
	vec devianceResiduals;

	ModelFit() {}

	// LEAST
	ModelFit( const bool & success, const vec &coef) : 
		success(success), coef(coef) {}

	// LOW
	ModelFit( const bool & success, const vec &coef, const vec &se, const double & dispersion, const double &rdf) : 
		success(success), coef(coef), se(se), dispersion(dispersion), rdf(rdf) {

		init();
	}

	// LOW with scalar coef and se
	ModelFit( const bool & success, const double &coef, const double &se, const double & dispersion, const double &rdf): 
		success(success), coef(vec(1, fill::value(coef))), se(vec(1, fill::value(se))), dispersion(dispersion), rdf(rdf) {}

	// MEDIUM
	ModelFit( const bool & success, const vec &coef, const vec &se, const double & dispersion, const double &rdf, const mat & vcov) :
		success(success), coef(coef), se(se), dispersion(dispersion), rdf(rdf), vcov(vcov) {

		init();
	}

	// HIGH
	ModelFit( const bool & success, const vec &coef, const vec &se, const double & dispersion, const double &rdf, const mat & vcov, const vec &residuals) : 
		success(success), coef(coef), se(se), dispersion(dispersion), rdf(rdf), vcov(vcov), residuals(residuals) {

		init();
	}

	// MOST
	ModelFit( const bool & success, const vec &coef, const vec &se, const double & dispersion, const double &rdf, const mat & vcov, const vec &residuals, const vec &hatvalues) : 
		success(success), coef(coef), se(se), dispersion(dispersion), rdf(rdf), vcov(vcov), residuals(residuals), hatvalues(hatvalues) {

		init();
	}

	// MAX
	void setDevResids( const vec &dr, const vec &y, const vec &mu){

    	// transform from residuals.glm
    	// d.res <- sqrt(pmax((object$family$dev.resids)(y, mu, 
        //     wts), 0))
        // ifelse(y > mu, d.res, -d.res)
    	devianceResiduals = sqrt(pmax(dr, 0));
    	uvec idx = find(y <= mu);
    	devianceResiduals.elem(idx) = -1.0*devianceResiduals.elem(idx);
	}	

	void setFittedValues( const vec &mu_in){
		mu = mu_in;
	}

	private:
	void init(){
		// if model didn't succeed, set values to nan
		if( ! success ){
			se.fill(datum::nan);
			vcov.fill(datum::nan);
		
			// if vcov is empty, set to be se.n_elem x se.n_elem
			// with elements nan
			if( vcov.n_elem == 0){
				vcov = mat(se.n_elem, se.n_elem);
				vcov.fill(datum::nan);
			}
		}
	}
};


typedef vector<ModelFit> ModelFitList;

class ModelFitLMM : public ModelFit {      
  public:  
  // additional information needed beyond ModelFit
  double logLik;
  vec weights, ru, y;
  double delta, sigSq_g, sigSq_e;
  int iter;

  ModelFitLMM(){}

   // LEAST
  ModelFitLMM(const bool & success, 
              const double &logLik,
              const vec &weights,
              const vec &ru,
              const vec &y,
              const double &delta,
              const double &sigSq_g,
              const double &sigSq_e,
              const int &iter, 
              const vec &coef) : 
    ModelFit( success, coef),
    logLik(logLik), weights(weights), ru(ru), y(y), delta(delta), sigSq_g(sigSq_g), sigSq_e(sigSq_e), iter(iter)
    {} 

  // LOW
  ModelFitLMM(const bool & success, 
              const double &logLik,
              const vec &weights,
              const vec &ru,
              const vec &y,
              const double &delta,
              const double &sigSq_g,
              const double &sigSq_e,
              const int &iter, 
              const vec &coef, 
              const vec &se, 
              const double &rdf) : 
    ModelFit( success, coef, se, sigSq_e, rdf),
    logLik(logLik), weights(weights), ru(ru), y(y), delta(delta), sigSq_g(sigSq_g), sigSq_e(sigSq_e), iter(iter)
    {} 
  
  // MEDIUM
  ModelFitLMM(const bool & success, 
              const double &logLik,
              const vec &weights,
              const vec &ru,
              const vec &y,
              const double &delta,
              const double &sigSq_g,
              const double &sigSq_e,
              const int &iter, 
              const vec &coef, 
              const vec &se, 
              const double &rdf, 
              const mat & vcov) : 
    ModelFit( success, coef, se, sigSq_e, rdf, vcov),
    logLik(logLik), weights(weights), ru(ru), y(y), delta(delta), sigSq_g(sigSq_g), sigSq_e(sigSq_e), iter(iter) 
    {} 

  // HIGH
  ModelFitLMM(const bool & success, 
              const double &logLik,
              const vec &weights,
              const vec &ru,
              const vec &y,
              const double &delta,
              const double &sigSq_g,
              const double &sigSq_e,
              const int &iter, 
              const vec &coef, 
              const vec &se, 
              const double &rdf, 
              const mat & vcov, 
              const vec &residuals) : 
    ModelFit( success, coef, se, sigSq_e, rdf, vcov, residuals),
    logLik(logLik), weights(weights), ru(ru), y(y), delta(delta), sigSq_g(sigSq_g), sigSq_e(sigSq_e), iter(iter)
    {}   

  // MOST
  ModelFitLMM(const bool & success, 
              const double &logLik,
              const vec &weights,
              const vec &ru,
              const vec &y,
              const double &delta,
              const double &sigSq_g,
              const double &sigSq_e,
              const int &iter, 
              const vec &coef, 
              const vec &se, 
              const double &rdf, 
              const mat & vcov, 
              const vec &residuals, 
              const vec &hatvalues) : 
    ModelFit( success, coef, se, sigSq_e, rdf, vcov, residuals, hatvalues),
    logLik(logLik), weights(weights), ru(ru), y(y), delta(delta), sigSq_g(sigSq_g), sigSq_e(sigSq_e), iter(iter)
    {}   
};

typedef vector<ModelFitLMM> ModelFitLMMList;

}

#endif