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
// #include "glm_family.h"
#include "glm.h"
// #include "linearRegression.h"

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
            const ModelDetail md = LOW);

	// extract results
  ModelFitLMM get_result();

  private:
  fastlmm<T1,T2,T3> fit;

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
			            const ModelDetail md){

	Rcpp::Rcout << "fastglmm..." << std::endl;

	// Initialize eta
	ModelFitGLM fit_init = GLM(X, y, family, md, weights, offset, nullptr, {}, 1e-4, 5);
	ModelFit fit_init2 = lm(X, y, md);

	// vec etc = fit_init->eta;

	// PQL iterations
	for(int i=0; i<10; i++){

		// update mu, zz, wz, eta, 

		// fit fastlmm
		fit = fastlmm(y, X, U, s, weights, MAX);
	}


	

}



template <typename T1, typename T2, typename T3> 
ModelFitLMM fastglmm<T1, T2, T3>::get_result(){

	return fit.get_result();




}

} // end namespace
#endif