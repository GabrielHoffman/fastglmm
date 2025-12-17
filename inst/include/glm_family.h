/***************************************************************
 * @file	glm_family.h
 * @author	Gabriel Hoffman
 * @email	gabriel.hoffman@mssm.edu
 * @brief	Define GLM link functions
 * Copyright (C) 2024 Gabriel Hoffman
 **************************************************************/


#ifndef _GLM_FAMILY_H_
#define _GLM_FAMILY_H_

#include <string>
#include <regex>
#include "misc.h"

using namespace arma;	
using namespace std;


namespace fastglmmLib {

/** base class
 * This should be a virtual class
 * Leaving these empty or using override in the derived
 * classes causes:
 * 
 * symbol not found in flat namespace '__ZN9GLMFamilyD2E
*/
class GLMFamily {
	public: 
	GLMFamily() {}

	virtual ~GLMFamily() {};
	virtual vec link( const vec &mu) const {return vec(1);}
	virtual vec linkinv( const vec &eta) const {return vec(1);}
	virtual vec mu_eta( const vec &eta) const {return vec(1);}
	virtual vec variance( const vec &mu) const {return vec(1);}

	// compute deviance redisuals as in gaussian()$dev.resids
	// But also need to transform as in residuals.glm(): 
	// in ModelFit::setDevResids()
	virtual vec dev_resids( const vec &y, const vec &mu, const vec &weights) const {return vec(1);}
	virtual vec initialize( const vec &y, const vec &weights) const {return vec(1);}
	virtual bool estimateDispersion() const {return true;}
	virtual string family() const {return "GLMFamily";};	
	virtual bool isCountModel() const {return true;}
	virtual void setOverdispersion( const double & value){}
};

class GaussianIdentity :
	virtual public GLMFamily {

	public:
	GaussianIdentity() {}

	~GaussianIdentity() {}

	vec link( const vec &mu) const {
		return mu;
	}
	vec linkinv( const vec &eta) const {
		return eta;
	}
	vec mu_eta( const vec &eta) const {
		return vec(eta.n_elem, fill::ones);
	}
	vec variance( const vec &mu) const {
		return vec(mu.n_elem, fill::ones);
	}
	vec dev_resids( const vec &y, const vec &mu, const vec &weights) const {
		// wt * ((y - mu)^2)
		return weights % square(y-mu);	
	}
	vec initialize( const vec &y, const vec &weights) const {
		return y;
	}
	bool estimateDispersion() const {return true;}
	string family() const {return "GaussianIdentity";}
	bool isCountModel() const { return false; }
};

class BinomialLogit :
	virtual public GLMFamily {

	public:
	BinomialLogit() {}

	~BinomialLogit() {}

	vec link( const vec &mu) const {
		return log(mu / (1-mu));
	}
	vec linkinv( const vec &eta) const {
		return 1.0 / (1.0 + exp(-1.0*eta));
	}
	vec mu_eta( const vec &eta) const {
		vec v = exp(-1.0*eta);
		return v / square( 1.0 + v);
	}
	vec variance( const vec &mu) const {
		return mu % (1.0 - mu);
	}	
	vec dev_resids( const vec &y, const vec &mu, const vec &weights) const {

		// 2 * rwt[i] * (y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
		return 2.0 * weights % (y_log_y(y, mu) + y_log_y(1.0 - y, 1.0 - mu));
	}
	vec initialize( const vec &y, const vec &weights) const {
		return (weights % y + 0.5)/(weights + 1.0);
	}
	bool estimateDispersion() const {return false;}
	string family() const {return "BinomialLogit";}
	bool isCountModel() const { return false; }
};


class QuasibinomialLogit :
	public BinomialLogit {
	bool estimateDispersion() const {return true;}
	string family() const {return "QuasibinomialLogit";}
};


class BinomialProbit :
	virtual public GLMFamily {

	public:
	BinomialProbit(){}

	~BinomialProbit() {}

	vec link( const vec &mu) const {
		// qnorm(mu)
		return qnorm(mu);
	}
	vec linkinv( const vec &eta) const { 
	  // pnorm(eta)
	  return normcdf( pmin(pmax(eta, -thresh), thresh) );
	}
	vec mu_eta( const vec &eta) const {
		// pmax(dnorm(eta), .Machine$double.eps)
		return pmax(normpdf(eta), tol);
	}
	vec variance( const vec &mu) const {
		return mu % (1.0 - mu);
	}
	vec dev_resids( const vec &y, const vec &mu, const vec &weights) const {
		
		// 2 * rwt[i] * (y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
		return 2.0 * weights % (y_log_y(y, mu) + y_log_y(1 - y, 1 - mu));
	}
	vec initialize( const vec &y, const vec &weights) const {
		return (weights % y + 0.5)/(weights + 1.0);
	}
	bool estimateDispersion() const {return false;}
	string family() const {return "BinomialProbit";}
	bool isCountModel() const { return false; }

	private:
	double tol = 2.220446e-16;
	// double thresh = -1.0 * R::qnorm( tol, 0.0, 1.0, 1, 0 );
	double thresh = 8.125891; // -1*qnorm(.Machine$double.eps)
};


class PoissonLog :
	virtual public GLMFamily {

	public:
	PoissonLog(){}

	~PoissonLog() {}

	vec link( const vec &mu) const {
		return log(mu);
	}
	vec linkinv( const vec &eta) const { 
	  // pmax(exp(eta), .Machine$double.eps)
	  return pmax(exp(eta), tol);
	}
	vec mu_eta( const vec &eta) const {
	  // pmax(exp(eta), .Machine$double.eps)
	  return pmax(exp(eta), tol);
	}
	vec variance( const vec &mu) const {
		return mu;
	}
	vec dev_resids( const vec &y, const vec &mu, const vec &weights) const {
		// r <- mu * wt
	  // p <- which(y > 0)
	  // r[p] <- (wt * (y * log(y/mu) - (y - mu)))[p]
	  // 2 * r

		vec res = mu % weights;
		uvec idx = find(y > 0.0);
		vec tmp = (weights % (y % log(y/mu) - (y - mu)));
		res.elem(idx) = tmp.elem(idx);

		return 2.0 * res;
	}
	vec initialize( const vec &y, const vec &weights) const {
		return y + 0.1;
	}
	bool estimateDispersion() const {return false;}
	string family() const {return "PoissonLog";}
	bool isCountModel() const { return true; }

	private:
	double tol = 2.220446e-16;
};


class QuasipoissonLog :
	public PoissonLog {
	bool estimateDispersion() const {return true;}
};

class NB :
	virtual public GLMFamily {

	public:	
	NB() {}

	NB(const double &theta) : 
		theta(theta) {}

	~NB() {}

	vec link( const vec &mu) const {
		return log(mu);
	}
	vec linkinv( const vec &eta) const { 
	  // pmax(exp(eta), .Machine$double.eps)
	  return pmax(exp(eta), tol);
	}
	vec mu_eta( const vec &eta) const {
	  // pmax(exp(eta), .Machine$double.eps)
	  return pmax(exp(eta), tol);
	}
	vec variance( const vec &mu) const {
		return mu + square(mu) / theta;
	}
	vec dev_resids( const vec &y, const vec &mu, const vec &weights) const {
		// 2 * wt * (y * log(pmax(1, y)/mu) - (y + .Theta) * log((y + .Theta)/(mu + .Theta)))

		return 2.0 * weights % (y % log(pmax(y, 1.0)/mu) - (y + theta) % log((y + theta)/(mu + theta)));
	}
	vec initialize( const vec &y, const vec &weights) const {
		// y + (y == 0)/6
		return y + accu(y == 0) / 6.0;
	}
	bool estimateDispersion() const {return true;}
	string family() const {return "NB";}
	bool isCountModel() const { return true; }

	void setOverdispersion( const double &value){
		theta = value;
	}

	double theta = std::numeric_limits<double>::quiet_NaN();

	private:
	double tol = 2.220446e-16;
};



/** Get smart pointer to GLMFamily object
 * @param family type of GLM: "gaussian", "gaussian/identity", "binomial/logit", "binomial/probit", "poisson/log", "quasibinomial/logit", "quasipoisson/log", "nb", nb:x" where x is a numeric value of theta  
 */ 
static shared_ptr<GLMFamily> getGLMFamily( const string &family){
	
	shared_ptr<GLMFamily> fam;

	if( family == "gaussian" || family == "gaussian/identity" ){
		fam = shared_ptr<GLMFamily>(new GaussianIdentity());
	}
	else if( family == "binomial/logit" ){
		fam = shared_ptr<GLMFamily>(new BinomialLogit());
	}
	else if( family == "binomial/probit" ){
		fam = shared_ptr<GLMFamily>(new BinomialProbit());
	}
	else if( family == "poisson/log" ){
		fam = shared_ptr<GLMFamily>(new PoissonLog());
	}
	else if( family == "quasibinomial/logit" ){
		fam = shared_ptr<GLMFamily>(new QuasibinomialLogit());
	}
	else if( family == "quasipoisson/log" ){
		fam = shared_ptr<GLMFamily>(new QuasipoissonLog());
	}
	else if( family == "nb" ){
		fam = shared_ptr<GLMFamily>(new NB());
	}
	else if( regex_search( family, regex("^nb:")) ){
		string theta_str = regex_replace( family, regex("^nb:"), "");
		double theta = atof(theta_str.c_str());
		fam = shared_ptr<GLMFamily>(new NB(theta));
	}else{
		throw invalid_argument( "Invalid GLM family: " + family );
	}

	return fam;
}

/** Return true of family is Poisson, quasi-Poisson or NB model
 */ 
static bool isCountModel(const string &family){

	shared_ptr<GLMFamily> fam = getGLMFamily( family );
	string famStr = fam->family();

	bool value = false;
	if( famStr == "PoissonLog" ||
		famStr == "QuasipoissonLog" || 
		famStr == "NB"){
		value = true;
	}

	return value;
}

}

#endif