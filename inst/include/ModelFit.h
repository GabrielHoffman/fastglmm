/***************************************************************
 * @file		ModelFit.h
 * @author	Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		Store parameters from model fit
 * Copyright (C) 2024 Gabriel Hoffman
 **************************************************************/

#ifndef MODEL_FIT_H_
#define MODEL_FIT_H_

#include <vector>
#include <string>
#include <regex>

#include "misc.h"

#include "glm_family.h"

using namespace arma;
using namespace std;

namespace fastglmmLib {

// Specify level of model detail to return from regression fit
typedef enum {
	LEAST,  // just beta
  LOW, 	  // baseline parameters: beta, se, dispersion, rdf
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
	double dispersion = datum::nan;
	double rdf = datum::nan;
	string ID;
	mat vcov;
	vec residuals;
	vec hatvalues;
	vec mu;
	vec devianceResiduals;
  double varFitted = datum::nan;

	ModelFit() {}

	// LEAST
	ModelFit( const bool & success, 
            const vec &coef) : 
		success(success), coef(coef) {}

	// LOW
	ModelFit( const bool & success, 
            const vec &coef, 
            const vec &se, 
            const double & dispersion, 
            const double &rdf) : 
		success(success), coef(coef), se(se), dispersion(dispersion), rdf(rdf) {

		init();
	}

	// LOW with scalar coef and se
	ModelFit( const bool & success, 
            const double &coef, 
            const double &se, 
            const double & dispersion, 
            const double &rdf): 
		success(success), coef(vec(1, fill::value(coef))), se(vec(1, fill::value(se))), dispersion(dispersion), rdf(rdf) {}

	// MEDIUM
	ModelFit( const bool & success, 
            const vec &coef, 
            const vec &se, 
            const double & dispersion, 
            const double &rdf, 
            const mat & vcov) :
		success(success), coef(coef), se(se), dispersion(dispersion), rdf(rdf), vcov(vcov) {

		init();
	}

	// HIGH
	ModelFit( const bool & success, 
            const vec &coef, 
            const vec &se, 
            const double & dispersion, 
            const double &rdf, 
            const mat & vcov, 
            const vec &residuals) : 
		success(success), coef(coef), se(se), dispersion(dispersion), rdf(rdf), vcov(vcov), residuals(residuals) {

		init();
	}

	// MOST
	ModelFit( const bool & success, 
            const vec &coef, 
            const vec &se, const 
            double & dispersion, 
            const double &rdf, 
            const mat & vcov, 
            const vec &residuals, 
            const vec &hatvalues) : 
		success(success), coef(coef), se(se), dispersion(dispersion), rdf(rdf), vcov(vcov), residuals(residuals), hatvalues(hatvalues) {

		init();
	}

	// MAX
	void setDevResids( 
    const vec &dr, 
    const vec &y, 
    const vec &mu,
    const vec w = {}){

  	// transform from residuals.glm
  	// d.res <- sqrt(pmax((object$family$dev.resids)(y, mu, 
    //     wts), 0))
    // ifelse(y > mu, d.res, -d.res)
  	devianceResiduals = sqrt(pmax(dr, 0));
  	uvec idx = find(y <= mu);
  	devianceResiduals.elem(idx) = -1.0*devianceResiduals.elem(idx);

    if( ! w.is_empty() ){
      // if weight is zero, set element to NAN
      devianceResiduals.elem(find(w == 0.0)).fill(datum::nan);
    }
	}	

	void setFittedValues( const vec &mu_in, const vec w = {}){
		mu = mu_in;

    if( ! w.is_empty() ){
      // if weight is zero, set element to NAN
      mu.elem(find(w == 0.0)).fill(datum::nan);
    }
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


class ModelFitGLM : public ModelFit {

  public:

  ModelFitGLM() {}

  ModelFitGLM( ModelFit &gmf, const string &family, const int &niter ) :
    ModelFit(gmf), family(family), niter(niter) {

    // extract theta values from "nb:theta"
    if( regex_search( family, regex("^nb:")) ){
      string theta_str = regex_replace( family, regex("^nb:"), "");
      theta = atof(theta_str.c_str());
    }else{
      theta = datum::nan;
    }

    // Account for loss of degrees of freedom when predicted counts are zero in a count model
      // https://doi.org/10.1515/sagmb-2017-0010
      // if family is Poisson, quasipoisson or nb
    // Compute number of entries were predicted values are 
    // effectively zero counts
    shared_ptr<GLMFamily> fam = getGLMFamily( family );
    if( fam->isCountModel() ){
      nZeroPrediction = sum(gmf.mu < 1e-4);
    }
  }

  string family = "";
  double theta = datum::nan;
  double mu_mean = datum::nan;
  double y_mean = datum::nan;
  int niter = 0;
  double nZeroPrediction = 0;
};

class ModelFitLMM : public ModelFit {      
  public:  
  // additional information needed beyond ModelFit
  double logLik;
  vec weights, ru, y;
  double delta, sigSq_g, sigSq_e;
  int iter;  
  double w_mean; // Store mean of weights before scaling

  bool isSet_U = false;
  bool isSet_Usp = false;
  bool isSet_V = false;
  bool isSet_Vsp = false;
  mat U, V;
  sp_mat Usp, Vsp; 
  vec s;  

  // Precomputed values for Sattherthwaite DDF
  mat A_sat, B_sat, hessian_vc;

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
              const double &w_mean,
              const vec &coef) : 
    ModelFit( success, coef),
    logLik(logLik), weights(weights), ru(ru), y(y), delta(delta), sigSq_g(sigSq_g), sigSq_e(sigSq_e), iter(iter), w_mean(w_mean)
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
              const double &w_mean,
              const vec &coef, 
              const vec &se, 
              const double &rdf) : 
    ModelFit( success, coef, se, sigSq_e, rdf),
    logLik(logLik), weights(weights), ru(ru), y(y), delta(delta), sigSq_g(sigSq_g), sigSq_e(sigSq_e), iter(iter), w_mean(w_mean)
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
              const double &w_mean,
              const vec &coef, 
              const vec &se, 
              const double &rdf, 
              const mat & vcov) : 
    ModelFit( success, coef, se, sigSq_e, rdf, vcov),
    logLik(logLik), weights(weights), ru(ru), y(y), delta(delta), sigSq_g(sigSq_g), sigSq_e(sigSq_e), iter(iter), w_mean(w_mean) 
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
              const double &w_mean,
              const vec &coef, 
              const vec &se, 
              const double &rdf, 
              const mat & vcov, 
              const vec &residuals) : 
    ModelFit( success, coef, se, sigSq_e, rdf, vcov, residuals),
    logLik(logLik), weights(weights), ru(ru), y(y), delta(delta), sigSq_g(sigSq_g), sigSq_e(sigSq_e), iter(iter), w_mean(w_mean)
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
              const double &w_mean,
              const vec &coef, 
              const vec &se, 
              const double &rdf, 
              const mat & vcov, 
              const vec &residuals, 
              const vec &hatvalues) : 
    ModelFit( success, coef, se, sigSq_e, rdf, vcov, residuals, hatvalues),
    logLik(logLik), weights(weights), ru(ru), y(y), delta(delta), sigSq_g(sigSq_g), sigSq_e(sigSq_e), iter(iter), w_mean(w_mean)
    {}   

    void setUS(const mat &U_, const vec &s_, const mat &V_){
      U = U_;
      V = V_;
      s = s_;
      isSet_U = true;
      isSet_V = true;
    }

    void setUS(const sp_mat &U_, const vec &s_, const mat &V_){
      Usp = U_;
      V = V_;
      s = s_;
      isSet_Usp = true;
      isSet_V = true;
    }

    void setUS(const mat &U_, const vec &s_, const sp_mat &V_){
      U = U_;
      Vsp = V_;
      s = s_;
      isSet_U = true;
      isSet_Vsp = true;
    }

    void setUS(const sp_mat &U_, const vec &s_, const sp_mat &V_){
      Usp = U_;
      Vsp = V_;
      s = s_;
      isSet_Usp = true;
      isSet_Vsp = true;
    }

    void set_w_mean( const double &value){
      w_mean = value;
    }
};


class ModelFitGLMM : public ModelFitLMM {

  public:

  ModelFitGLMM() {}

  ModelFitGLMM( ModelFitLMM &gmf, const string &family, const int &iter ) :
    ModelFitLMM(gmf), family(family) {

    this->iter = iter;
    
    // extract theta values from "nb:theta"
    if( regex_search( family, regex("^nb:")) ){
      string theta_str = regex_replace( family, regex("^nb:"), "");
      theta = atof(theta_str.c_str());
    }else{
      theta = datum::nan;
    }

    // Account for loss of degrees of freedom when predicted counts are zero in a count model
      // https://doi.org/10.1515/sagmb-2017-0010
      // if family is Poisson, quasipoisson or nb
    // Compute number of entries were predicted values are 
    // effectively zero counts
    shared_ptr<GLMFamily> fam = getGLMFamily( family );
    if( fam->isCountModel() ){
      nZeroPrediction = sum(gmf.mu < 1e-4);
    } 
  }

  vec blup;
  string family = "";
  double theta = datum::nan;
  double mu_mean = datum::nan;
  double y_mean = datum::nan;
  double nZeroPrediction = 0;
};


typedef vector<ModelFit> ModelFitList;
typedef vector<ModelFitGLM> ModelFitGLMList;
typedef vector<ModelFitLMM> ModelFitLMMList;
typedef vector<ModelFitGLMM> ModelFitGLMMList;



}

#endif