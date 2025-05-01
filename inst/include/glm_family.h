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


/** Virtual base class
*/
class GLMFamily {
	public: 
	GLMFamily() {}
	virtual ~GLMFamily() = 0;
	virtual vec link( const vec &mu) const = 0;
	virtual vec linkinv( const vec &eta) const = 0;
	virtual vec mu_eta( const vec &eta) const = 0;
	virtual vec variance( const vec &mu) const = 0;

	// compute deviance redisuals as in gaussian()$dev.resids
	// But also need to transform as in residuals.glm(): 
	// in ModelFit::setDevResids()
	virtual vec dev_resids( const vec &y, const vec &mu, const vec &weights) const = 0;
	virtual vec initialize( const vec &y, const vec &weights) const = 0;	
	virtual bool estimateDispersion() const = 0;
	virtual string family() const = 0;	
	virtual bool isCountModel() const = 0;
};

namespace fastglmmLib {
class GaussianIdentity :
	public virtual GLMFamily {

	public:
	GaussianIdentity() {}

	~GaussianIdentity() override {}

	vec link( const vec &mu) const override {
		return mu;
	}
	vec linkinv( const vec &eta) const override {
		return eta;
	}
	vec mu_eta( const vec &eta) const override {
		return vec(eta.n_elem, fill::ones);
	}
	vec variance( const vec &mu) const override {
		return vec(mu.n_elem, fill::ones);
	}
	vec dev_resids( const vec &y, const vec &mu, const vec &weights) const override {
		// wt * ((y - mu)^2)
		return weights % pow(y-mu,2);	
	}
	vec initialize( const vec &y, const vec &weights) const override {
		return y;
	}
	bool estimateDispersion() const override {return true;}
	string family() const override {return "GaussianIdentity";}
	bool isCountModel() const override { return false; }
};

class BinomialLogit :
	public virtual GLMFamily {

	public:
	BinomialLogit() {}

	~BinomialLogit() override {}

	vec link( const vec &mu) const override {
		return log(mu / (1-mu));
	}
	vec linkinv( const vec &eta) const override {
		return 1.0 / (1.0 + exp(-1.0*eta));
	}
	vec mu_eta( const vec &eta) const override {
		vec v = exp(-1.0*eta);
		return v / pow( 1.0 + v, 2);
	}
	vec variance( const vec &mu) const override {
		return mu % (1.0 - mu);
	}	
	vec dev_resids( const vec &y, const vec &mu, const vec &weights) const override {

		// 2 * rwt[i] * (y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
		return 2.0 * weights % (y_log_y(y, mu) + y_log_y(1.0 - y, 1.0 - mu));
	}
	vec initialize( const vec &y, const vec &weights) const override {
		return (weights % y + 0.5)/(weights + 1.0);
	}
	bool estimateDispersion() const override {return false;}
	string family() const override {return "BinomialLogit";}
	bool isCountModel() const override { return false; }
};


class QuasibinomialLogit :
	public BinomialLogit {
	bool estimateDispersion() const override {return true;}
	string family() const override {return "QuasibinomialLogit";}
};


class BinomialProbit :
	public virtual GLMFamily {

	public:
	BinomialProbit(){}

	~BinomialProbit() override {}

	vec link( const vec &mu) const override {
		// qnorm(mu)
		return qnorm(mu);
	}
	vec linkinv( const vec &eta) const override { 
	    // pnorm(eta)
	    return normcdf( pmin(pmax(eta, -thresh), thresh) );
	}
	vec mu_eta( const vec &eta) const override {
		// pmax(dnorm(eta), .Machine$double.eps)
		return pmax(normpdf(eta), tol);
	}
	vec variance( const vec &mu) const override {
		return mu % (1.0 - mu);
	}
	vec dev_resids( const vec &y, const vec &mu, const vec &weights) const override {
		
		// 2 * rwt[i] * (y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
		return 2.0 * weights % (y_log_y(y, mu) + y_log_y(1 - y, 1 - mu));
	}
	vec initialize( const vec &y, const vec &weights) const override {
		return (weights % y + 0.5)/(weights + 1.0);
	}
	bool estimateDispersion() const override {return false;}
	string family() const override {return "BinomialProbit";}
	bool isCountModel() const override { return false; }

	private:
	double tol = 2.220446e-16;
	// double thresh = -1.0 * R::qnorm( tol, 0.0, 1.0, 1, 0 );
	double thresh = 8.125891; // -1*qnorm(.Machine$double.eps)
};


class PoissonLog :
	public virtual GLMFamily {

	public:
	PoissonLog(){}

	~PoissonLog() override {}

	vec link( const vec &mu) const override {
		return log(mu);
	}
	vec linkinv( const vec &eta) const override { 
	    // pmax(exp(eta), .Machine$double.eps)
	    return pmax(exp(eta), tol);
	}
	vec mu_eta( const vec &eta) const override {
	    // pmax(exp(eta), .Machine$double.eps)
	    return pmax(exp(eta), tol);
	}
	vec variance( const vec &mu) const override {
		return mu;
	}
	vec dev_resids( const vec &y, const vec &mu, const vec &weights) const override {
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
	vec initialize( const vec &y, const vec &weights) const override {
		return y + 0.1;
	}
	bool estimateDispersion() const override {return false;}
	string family() const override {return "PoissonLog";}
	bool isCountModel() const override { return true; }

	private:
	double tol = 2.220446e-16;
};


class QuasipoissonLog :
	public PoissonLog {
	bool estimateDispersion() const override {return true;}
};

class NB :
	public virtual GLMFamily {

	public:
	NB(const double &theta) : 
		theta(theta) {}

	~NB() override {}

	vec link( const vec &mu) const override {
		return log(mu);
	}
	vec linkinv( const vec &eta) const override { 
	    // pmax(exp(eta), .Machine$double.eps)
	    return pmax(exp(eta), tol);
	}
	vec mu_eta( const vec &eta) const override {
	    // pmax(exp(eta), .Machine$double.eps)
	    return pmax(exp(eta), tol);
	}
	vec variance( const vec &mu) const override {
		return mu + pow(mu, 2) / theta;
	}
	vec dev_resids( const vec &y, const vec &mu, const vec &weights) const override {
		// 2 * wt * (y * log(pmax(1, y)/mu) - (y + .Theta) * log((y + .Theta)/(mu + .Theta)))

		return 2.0 * weights % (y % log(pmax(y, 1.0)/mu) - (y + theta) % log((y + theta)/(mu + theta)));
	}
	vec initialize( const vec &y, const vec &weights) const override {
		// y + (y == 0)/6
		return y + accu(y == 0) / 6.0;
	}
	bool estimateDispersion() const override {return true;}
	string family() const override {return "NB";}
	bool isCountModel() const override { return true; }

	double theta;

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
	else if( regex_search( family, regex("^nb:")) ){
		string theta_str = regex_replace( family, regex("^nb:"), "");
		double theta = atof(theta_str.c_str());
		fam = shared_ptr<GLMFamily>(new NB(theta));
	}else{
		throw invalid_argument( "Invalid GLM family: " + family );
	}

	return fam;
}


}

#endif