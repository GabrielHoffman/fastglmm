/***************************************************************
 * @file		fastglmm_fit.h
 * @author	Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		Fit generalizd linear mixed model
 * Copyright (C) 2024 Gabriel Hoffman
 **************************************************************/

#ifndef _FASTGLMM_FIT_H_
#define _FASTGLMM_FIT_H_

// if -D USE_R, use RcppArmadillo library
#ifdef USE_R
// [[Rcpp::depends(RcppParallel, roptim)]]	
#include <RcppArmadillo.h>
#include <roptim.h>
#else
#include <armadillo>
#endif

#include "fastlmm_fit.h"
#include "glm_family.h"
#include "glm.h"
#include "ModelFit.h"
#include "spectralDecomp.h"
#include "GLMMLaplace.h"

using namespace arma;
using namespace std;
using namespace fastglmmLib;

namespace fastglmmLib {

typedef enum {
  PQL,
  LAPLACE
} GLMM_Method;

// Order of template variables
// T1 y
// T2 X
// T3 U
template <typename T1, typename T2, typename T3> 
class fastglmm {       
  public:  

	// constructor, minimal
	fastglmm(){};

	fastglmm(	const T1 &y, 
						const T2 &X, 
						const spectralDecomp<T3> &dcmp,
						const vec &weights,
						const vec &offset,
						const string &family, 
						const ModelDetail md = LOW, 
						const GLMM_Method &method = PQL,
						const double &tol = 1e-5,
						const double &tol_eta = 1e-7,
						const int &maxit = 100,
						const double &lambda = 0,
						const double &delta = -1,
						const double &left = -10,
						const double &right = 10,
						const bool &returnUS = false,
						const bool &doCoxReid = false);

	const vec residuals() const ; // Response
	const vec residuals_pearson() const ; // Pearson
	const vec fitted() const ;
	const vec devianceResiduals() const;
	const vec get_blup() const;
	const GLMM_Method get_method() const {
		return method;
	};

	// extract results
  ModelFitGLMM get_result();

  private:
  void _fastglmm_pql(const bool &estimateTheta, GLMWork *work);
  void _fastglmm_laplace(vec &betaInit, const bool &estimateTheta, GLMWork *work);

  fastlmm<T1,T2,T3> fit;
  vec y;
  mat X;
  vec blup;
  vec weights, mu, eta, offset; 
	spectralDecomp<T3> dcmp;
  string family;
  double lambda;
  bool returnUS;
  int niter_pql;
  double w_mean; 
  double mu_mean = datum::nan;
  double y_mean = datum::nan;
  double eta_var = datum::nan;
  ModelDetail md;
	GLMM_Method method;
	shared_ptr<GLMFamily> fam;
	bool isValid = true;
	mat V;
	double theta, delta;
	double tol, tol_eta;
	double left, right;
	int maxit;
	bool doCoxReid;
  unique_ptr<ModelFitLMM> mf_laplace;
};

template <typename T1, typename T2, typename T3> 
fastglmm<T1, T2, T3>::fastglmm(
	const T1 &y, 
	const T2 &X, 
	const spectralDecomp<T3> &dcmp,
	const vec &weights,
	const vec &offset,
	const string &family, 
	const ModelDetail md, 
	const GLMM_Method &method,
	const double &tol,
	const double &tol_eta,
	const int &maxit,
	const double &lambda,
	const double &delta,
	const double &left,
	const double &right,
	const bool &returnUS,
	const bool &doCoxReid):
	y(y), 
	X(X),
	weights(weights), 
	offset(offset), 
	dcmp(dcmp), 
	family(family), 
	lambda(lambda),
	returnUS(returnUS),
	md(md),
	method(method),
	delta(delta),
	tol(tol),
	tol_eta(tol_eta),
	left(left),
	right(right),
	maxit(maxit),
	doCoxReid(doCoxReid)
	{

	fam = getGLMFamily( family );

	// if Negative Binomial with unspecified theta
	// estimate theta, and initialize with Poisson GLM
	bool estimateTheta = family == "nb" ? true : false;
	if( estimateTheta ){
		this->family = "poisson/log";
	}

	checkResponse(y, this->family);

	GLMWork *work = new GLMWork();

	// Initialize eta
	// just need a rough starting value
	ModelFitGLM fit_init = GLM(X, y, this->family, LEAST, weights, offset, work, {}, 1e-2, 3, lambda);

	// if GLM fails, 
	if( ! fit_init.success ){

		// fit simple fastlmm, 
		// set maxit to zero to avoid PQL iterations
		// set isValid to FALSE
		work->w = vec(y.n_elem, fill::value(1));
			this->dcmp.reweight(work->w);
		fit = fastlmm(y, X, this->dcmp, work->w, LEAST, lambda);
		fit.eval_delta( 1 ); 
		this->maxit = 0;
		isValid = false;
	}

	// Estimation
  switch( method ){
    case PQL:
			_fastglmm_pql(estimateTheta, work);
			break;
		case LAPLACE:
			_fastglmm_laplace(fit_init.coef, estimateTheta, work);
			break;
	}

	if( estimateTheta ){
		// update family to include estimated theta
		this->family = "nb:" + to_string(theta);
	}

	// Use result of lmm and inverse link
	// to get final value of mu
	if( isValid ){
		if( method == PQL){
			eta = fit.fitted() + offset;
  		mu = fam->linkinv(eta);

			// save the BLUP from the fastlmm into the fastglmm object
			// apply V to correct ordering
			blup = this->dcmp.get_V() * fit.blup();
  	}
	}else{
		eta = vec(offset.n_elem, fill::value(datum::nan));
		mu = vec(offset.n_elem, fill::value(datum::nan));
		blup = vec(this->dcmp.get_V().n_rows, fill::value(datum::nan));
	}

	// Compute mean of mu
	// use robust mean to avoid influence of outliers
	// This can happen with many zero and a few large values
	mu_mean = robust_mean(mu, 4);

	eta_var = var(eta);

	delete work;
}


template <typename T1, typename T2, typename T3> 
const vec fastglmm<T1, T2, T3>::residuals() const {

	// Response residuals
	return( y - mu ); 
}


template <typename T1, typename T2, typename T3> 
const vec fastglmm<T1, T2, T3>::residuals_pearson() const {
	// Pearson residuals
	// (y - mu) * sqrt(wts) / sqrt(fam$variance(mu))
	return (y - mu) % sqrt(weights) / sqrt(fam->variance(mu));
}


template <typename T1, typename T2, typename T3> 
const vec fastglmm<T1, T2, T3>::fitted() const {

	vec v;

  switch( method ){
    case PQL:
			v = fam->linkinv( fit.fitted() );
			break;
		case LAPLACE:
			v = mu;
			break;
	}

	return v;
}


template <typename T1, typename T2, typename T3> 
const vec fastglmm<T1, T2, T3>::devianceResiduals() const {

	// transform from residuals.glm
	// d.res <- sqrt(pmax((object$family$dev.resids)(y, mu, 
  //     wts), 0))
  // ifelse(y > mu, d.res, -d.res)

	// compute raw deviance residuals
	vec dr = fam->dev_resids(y, mu, weights);

	vec drMod = sqrt(pmax(dr, 0));
	uvec idx = find(y <= mu);
	drMod.elem(idx) = -1.0*drMod.elem(idx);

	return drMod;
}


template <typename T1, typename T2, typename T3> 
ModelFitGLMM fastglmm<T1, T2, T3>::get_result(){

	ModelFitLMM res1;

	if( method == PQL){
		res1 = fit.get_result(returnUS);
		res1.set_w_mean( w_mean );
	}else{
		res1 = *mf_laplace;
		niter_pql = res1.iter;
	}

	// if model is not valid, it failed before final calculations
	// so set values to NAN matching ModelDetail
	if( ! isValid ){
		// number of coefs
		int p = X.n_cols;
		res1.coef = vec(p, fill::value(datum::nan));

		switch( md ){
	    case MAX:       
	    case MOST:
	      res1.mu.fill(datum::nan); 
	      res1.hatvalues.fill(datum::nan); 
	    case HIGH: 
	      res1.residuals.fill(datum::nan);
	    case MEDIUM: 
				res1.vcov = mat(p, p, fill::value(datum::nan));
	    case LOW: 
				res1.se = vec(p, fill::value(datum::nan));
				res1.rdf = datum::nan;
	    case LEAST: 
	      break;
		}
	}

	ModelFitGLMM mf(res1, family, niter_pql);

	mf.mu_mean = mu_mean;
	mf.y_mean = mean(y);
	mf.varFitted = eta_var;

	// BLUP
	mf.blup = blup;

  // QL dispersion based on Pearson residuals
  // sum((w*r^2)[w > 0]) / df.r
  // need to drop elements without nan  
  vec rp = residuals_pearson();
  vec rp_clean = rp.elem(find_finite(rp));
  double disp = dot(rp_clean, rp_clean) / res1.rdf ;

	if( fam->estimateDispersion() ){ 
	  mf.dispersion = disp;
	}else{
		mf.dispersion = 1.0;
		if( method == PQL ){
			// unscale variances by dispersion
			mf.vcov /= disp;
			mf.se /= sqrt(disp);
		}
	}

  if( md == MAX ){
		mf.devianceResiduals = devianceResiduals();
  }

  if( md >= HIGH ){
  	// response residuals
		mf.residuals = residuals(); 
  }
  
  if( md >= MOST ){
  	// Respones residuals
		mf.mu = fitted(); 
		// mf.hatvalues = hatvalues(); 
  } 

	return mf;
}


template <typename T1, typename T2, typename T3> 
void fastglmm<T1, T2, T3>::_fastglmm_pql(const bool &estimateTheta, GLMWork *work){

	vec eta_old;
	int iter_in = 0;
	uvec idx_drop = find(weights < 1e-12);
	double n_active = weights.n_elem - idx_drop.n_elem;

	CountTable ct;
	if( estimateTheta && maxit > 0){
		// Precompute lgamma() on each unique count
		ct = CreateLUT(y, weights);	
	}

	// PQL iterations
	for(niter_pql=0; niter_pql<maxit; niter_pql++){

		// if Negative Binomial with unspecified theta
		if( estimateTheta ){
			theta = nb_theta_ml(y, work->mu, y.n_elem, weights, X, doCoxReid, ct, -5, 20);
			fam->setOverdispersion( theta );
		}

		// update mu, eta, z, w, eta, 
		if( niter_pql == 0){
			work->eta = work->eta + offset;
		}else{
			eta_old = work->eta;
			work->eta = fit.fitted() + offset;

			// convergence criteria based on norm of eta change
			if( norm(work->eta - eta_old) < tol_eta){
				break;
			}
		}

		// entries with zero weights have NAN value
		work->eta.elem( idx_drop ).zeros();

		// mu <- family$linkinv(eta)
		work->mu = fam->linkinv( work->eta );

		// mu.eta.val <- family$mu.eta(eta)
		work->gprime = fam->mu_eta( work->eta );

		// zz <- eta + (y.orig - mu)/mu.eta.val - offset
		work->z = (work->eta - offset) + (y - work->mu) / work->gprime;

 		// wz <- w * mu.eta.val^2/family$variance(mu)
		if( fam->family() == "NB"){
			// Numerically stable evaluation for NB		
			theta = fam->getOverdispersion();		
	    vec x = work->eta - log(theta);
	    work->w = weights % (theta * (0.5 * (1.0 + tanh(0.5 * x))));	    
		}else{
			work->w = square(work->gprime) % (weights / fam->variance( work->mu ));
		}

		// wz <- wz / mean(wz)
		w_mean = sum(work->w) / n_active;
		work->w = work->w / w_mean;

		// if model has nan values in z or w, 
		// it can't be fit
		// so set beta values to nan 
		// and set isValid to false
		if( work->z.has_nan() || work->w.has_nan() ){
			fit.set_model_failure();
			isValid = false;
			break;
		}

		// recompute U and s since work->w changed
		this->dcmp.reweight(work->w);

		// fit fastlmm
		fit = fastlmm(work->z, X, this->dcmp, work->w, LEAST, lambda);

		if( delta > 0 ){
			fit.eval_delta( delta ); 
		}else{
			fit.estimate_delta(left, right, tol);	
		}

		// check condition number of vcov matrix
		if( fit.get_vcov_kappa() > 1e10 ){
			fit.set_model_failure();
			isValid = false;
			break;
		}

		// increment interation count
		iter_in += fit.get_iter();
	}
	

	// Final fit with ModelDetail md
	if( isValid && (md > LEAST) ){		

		// fit fastlmm
		fit = fastlmm(work->z, X, this->dcmp, work->w, md, lambda);

		if( delta > 0 ){
			fit.eval_delta( delta ); 
		}else{
			fit.estimate_delta(left, right, tol);	
		}
	}
}


template <typename T1, typename T2, typename T3> 
void fastglmm<T1, T2, T3>::_fastglmm_laplace(vec &betaInit, const bool &estimateTheta, GLMWork *work){

	GLMMLaplace<T1, T2, T3> obj(y, X, dcmp, weights, offset, lambda, "nb", LOW, work);

	if( ! estimateTheta ){
		theta = fam->getOverdispersion();
		obj.set_theta(theta);
	}

	// Use Nelder-Mead optimizer
	// since BFGS evalutes at -Inf 
  Roptim<GLMMLaplace<T1, T2, T3> > opt("Nelder-Mead");
  // opt.control.trace = 2;
  opt.set_hessian(true);  
  opt.control.reltol = tol / 5000;
	opt.control.maxit  = 5000;

  // set initial values
  vec param(betaInit.n_elem+1+estimateTheta);
  param.subvec(0, betaInit.n_elem-1) = betaInit; // beta
  param[betaInit.n_elem] = log(.5); // log_sigSq_b
  if( estimateTheta ){
  	param[betaInit.n_elem+1] = 2; // log_theta
  }

  opt.minimize(obj, param);

  // opt.print();

  // Get results:
  // split parameter vectors in separate variables
  vec beta = opt.par().subvec(0, betaInit.n_elem-1);
  double sigSq_g = pow(exp(opt.par()[betaInit.n_elem]), 2);
  if( estimateTheta ){
		theta = exp(opt.par()[betaInit.n_elem+1]);
	}
	int niter = opt.fncount();

	// Hessian includes sigSq_g and theta
	// subset to just get hessian of coefficients
	span sp(0, betaInit.n_elem-1);
	mat vcov = inv(opt.hessian()(sp, sp));
  vec se = sqrt(vcov.diag());

  this->mu = work->mu;
  this->eta = work->eta;
  this->weights = work->w;

  // BLUP
	// a = eta - X %*% beta - offset
	// BLUP = solve(crossprod(Z)) %*% crossprod(Z,a)
	T3 Z = obj.get_Z();
	T3 cpZ = Z.t() * Z;
	vec a = vec(cpZ.diag());
  vec BLUP = (Z.t() * (eta - X * beta - offset)) / a;

  // working residuals
  vec r = work->z - X * beta - offset - Z * BLUP;
  double sigSq_e = dot(r,r) / r.n_elem;
  delta = sigSq_e/sigSq_g;

  double rdf = y.n_elem - X.n_cols;
  fam->setOverdispersion( theta );

  mf_laplace = make_unique<ModelFitLMM>(true, -1.0*opt.value(), work->w, vec(), y, delta, sigSq_g, sigSq_e, niter, mean(work->w), beta, se, rdf, vcov, residuals() );

  T3 U = dcmp.get_U();

  // if returnUS
  // return U and s 
  if( returnUS ){
    switch( dcmp.get_type() ){
      case GENERAL:
        mf_laplace->setUS(U, dcmp.get_s(), dcmp.get_V());
        break;
      case CATEGORICAL:
        // V is identity
        mf_laplace->setUS(U, dcmp.get_s(), eye<sp_mat>(U.n_cols, U.n_cols));
        break;
    };
  }

  // // linearized response
  // work->eta = X * beta + offset;
  // work->mu = exp(work->eta);
  // work->z = work->eta + (y - work->mu) / work->mu;

	// mat Yu = dcmp.get_U().t() * work->z;
  mat Xu = dcmp.get_U().t() * X;
  mat Gamma_XX = X.t() * X - Xu.t() * Xu;
  vec inv_s_delta = 1 / (dcmp.get_s()+delta);
  mat inv_s_delta_Xu = scaleEachCol(Xu, inv_s_delta);

  // mf_laplace->ru = Yu - Xu * beta - dcmp.get_U().t() * offset;

  // mf_laplace->ru.brief_print("ru:");

  // with PQL, BLUP is 
  // (sqrt(s) % ru) / (s + delta_hat)
  // we computed BLUP independently
  // BLUP = (sqrt(s) % ru) / (s + delta_hat)
  // so delta = ((sqrt(s) * ru) - BLUP*s ) / BLUP;
  // vec delta_hat = ((sqrt(dcmp.get_s()) % mf_laplace->ru) - BLUP % dcmp.get_s()) / BLUP;

  // delta_hat.brief_print("delta_hat:");

	// Rcpp::Rcout << "delta_hat: " << mean(delta_hat) <<endl;

  // Precompute values for Satterthwaite DDF to be
  // used later with V and L specified
  Satterthwaite ddf_sat(y.n_elem, sigSq_g, sigSq_e, dcmp.get_s(), Xu, Gamma_XX, inv_s_delta, inv_s_delta_Xu);

  // save precomputed values
  mf_laplace->hessian_vc = ddf_sat.get_hessian();
  mf_laplace->A_sat = ddf_sat.get_A();
  mf_laplace->B_sat = ddf_sat.get_B();
}





} // end namespace
#endif