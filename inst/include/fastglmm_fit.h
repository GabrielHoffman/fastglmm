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
						const bool &returnUS = false);

	// extract results
  ModelFitGLMM get_result();

  private:
  fastlmm<T1,T2,T3> fit;
  string family;
  bool returnUS;
  int niter_pql;
  double w_mean; 

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
							const bool &returnUS):
							family(family), returnUS(returnUS) {

	shared_ptr<GLMFamily> fam = getGLMFamily( family );

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
	ModelFitGLM fit_init = GLM(X, y, this->family, md, weights, offset, work, {}, 1e-2, 3);

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
		fit = fastlmm(work->z, X, dcmp.get_vectors(), dcmp.get_values(), work->w, MAX);
		fit.estimate_delta(-10, 10, tol);	

		// TODO
		// 1) narrow delta search over time?
		iter_in += fit.get_iter();
	}

	if( estimateTheta ){
		// update family to include estimated theta
		this->family = "nb:" + to_string(theta);
	}

	delete work;
}



template <typename T1, typename T2, typename T3> 
ModelFitGLMM fastglmm<T1, T2, T3>::get_result(){

	ModelFitLMM res1 = fit.get_result(returnUS);
	res1.set_w_mean( w_mean );

	return ModelFitGLMM(res1, family, niter_pql);
}

} // end namespace
#endif