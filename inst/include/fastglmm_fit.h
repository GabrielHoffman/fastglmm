/***************************************************************
 * @file		fastglmm_fit.h
 * @author	Gabriel Hoffman
 * @email		gabriel.hoffman@mssm.edu
 * @brief		Fit generalizd linear mixed model
 * Copyright (C) 2024 Gabriel Hoffman
 **************************************************************/

#ifndef _FASTGLMM_FIT_H_
#define _FASTGLMM_FIT_H_

// if -D ARMA, use plain armadillo library
#ifdef ARMA
#include <armadillo>
#else
#include <RcppArmadillo.h>
#endif

#include "fastlmm_fit.h"
#include "glm_family.h"
#include "glm.h"
#include "ModelFit.h"
#include "spectralDecomp.h"

using namespace arma;
using namespace std;
using namespace fastglmmLib;

namespace fastglmmLib {

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
						const T3 &U, 
						const vec &s,
						const vec &weights,
						const vec &offset,
						const string &family, 
						const ModelDetail md = LOW, 
						const double &tol = 1e-4,
						const double &tol_eta = 1e-4,
						const int &maxit = 100,
						const double &delta = -1,
						const double &left = -10,
						const double &right = 10,
						const bool &returnUS = false);

	const vec residuals(); // Pearson
	const vec fitted();
	const vec devianceResiduals();

	// extract results
  ModelFitGLMM get_result();

  private:
  fastlmm<T1,T2,T3> fit;
  string family;
  bool returnUS;
  int niter_pql;
  double w_mean; 
  ModelDetail md;
  vec y, weights, mu; 
	shared_ptr<GLMFamily> fam;
};

template <typename T1, typename T2, typename T3> 
fastglmm<T1, T2, T3>::fastglmm(
							const T1 &y, 
							const T2 &X, 
							const T3 &U, 
							const vec &s,
							const vec &weights,
							const vec &offset,
							const string &family, 
							const ModelDetail md, 
							const double &tol,
							const double &tol_eta,
							const int &maxit,
							const double &delta,
							const double &left,
							const double &right,
							const bool &returnUS):
							y(y), weights(weights), family(family), returnUS(returnUS), md(md) {

	fam = getGLMFamily( family );

	// if Negative Binomial with unspecified theta
	// estimate theta, and initialize with Poisson GLM
	bool estimateTheta = family == "nb" ? true : false;
	if( estimateTheta ){
		this->family = "poisson/log";
	}

	checkResponse(y, this->family);

	GLMWork *work = new GLMWork();
	spectralDecomp<T3> dcmp;
	vec eta_old;
	T3 Z = scaleEachRow(U, sqrt(s));

	// Initialize eta
	// just need a rough starting value
	ModelFitGLM fit_init = GLM(X, y, this->family, LEAST, weights, offset, work, {}, 1e-2, 3);

	int iter_in = 0;
	double theta;

	// PQL iterations
	for(niter_pql=0; niter_pql<maxit; niter_pql++){

		// if Negative Binomial with unspecified theta
		if( estimateTheta ){
			theta = nb_theta_ml(y, work->mu, y.n_elem, weights, {}, false);
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

		// mu <- family$linkinv(eta)
		work->mu = fam->linkinv( work->eta );

		// mu.eta.val <- family$mu.eta(eta)
		work->gprime = fam->mu_eta( work->eta );

		// zz <- eta + (y.orig - mu)/mu.eta.val - offset
		work->z = (work->eta - offset) + (y - work->mu) / work->gprime;

		// wz <- w * mu.eta.val^2/family$variance(mu)
		work->w = pow(work->gprime,2) % (weights / fam->variance( work->mu ));

		// wz <- wz / mean(wz)
		w_mean = mean(work->w);
		work->w = work->w / w_mean;

		// recompute U and s since work->w changed
		dcmp.initWithIndicator(Z, work->w);

		// fit fastlmm
		fit = fastlmm(work->z, X, dcmp.get_vectors(), dcmp.get_values(), work->w, LEAST);

		if( delta > 0 ){
      fit.eval_delta( delta ); 
    }else{
			fit.estimate_delta(left, right, tol);	
		}

		// increment interation count
		iter_in += fit.get_iter();
	}

	// Final fit with ModelDetail md
	if( md > LEAST ){		
		// fit fastlmm
		fit = fastlmm(work->z, X, dcmp.get_vectors(), dcmp.get_values(), work->w, md);

		if( delta > 0 ){
      fit.eval_delta( delta ); 
    }else{
			fit.estimate_delta(left, right, tol);	
		}
	}

	if( estimateTheta ){
		// update family to include estimated theta
		this->family = "nb:" + to_string(theta);
	}

	// save value
  mu = work->mu;

	delete work;
}



template <typename T1, typename T2, typename T3> 
const vec fastglmm<T1, T2, T3>::residuals(){

	// (y - mu) * sqrt(wts)/sqrt(fam$variance(mu))
	return (y - mu) % sqrt(weights) / sqrt(fam->variance(mu));
}


template <typename T1, typename T2, typename T3> 
const vec fastglmm<T1, T2, T3>::fitted(){

	return fam->linkinv( fit.fitted() );
}



template <typename T1, typename T2, typename T3> 
const vec fastglmm<T1, T2, T3>::devianceResiduals(){

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

	ModelFitLMM res1 = fit.get_result(returnUS);
	res1.set_w_mean( w_mean );

	ModelFitGLMM mf(res1, family, niter_pql);

  if( md == MAX ){
		mf.devianceResiduals = devianceResiduals();
  }

  if( md >= HIGH ){
		mf.residuals = residuals();
  }

	return mf;
}


} // end namespace
#endif