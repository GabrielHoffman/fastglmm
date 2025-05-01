
#ifndef SPECTRAL_DECOMP_H_
#define SPECTRAL_DECOMP_H_

// if -D ARMA, use plain armadillo library
#ifdef ARMA
#include <armadillo>
#else
#include <RcppArmadillo.h>
#endif

#include "misc.h"

using namespace arma;

namespace fastglmmLib {

template <typename T> 
class spectralDecomp {

  public:
  	spectralDecomp(){};

  	// Initialize with indicator matrix Z;
    void initWithIndicator( const T &Z, const vec &weights);

    // Initialize with U and s from eigen decomp
    void initWithEigenDecomp( const T &U, const vec &s);

    // Initialize with U and s from eigen decomp, and weights
    void initWithEigenDecomp( const T &U, const vec &s, const vec &weights);

    T get_vectors(){return U;}
    vec get_values(){return s;}

  private:
   T U;
   vec s;
};

// Initialize with indicator matrix
template <typename T> 
void spectralDecomp<T>::initWithIndicator( const T &Z, const vec &weights){

	s = (weights.t() * Z).t();
	U = scaleRowsCols(Z, sqrt(weights), 1 / sqrt(s));
}


// Initialize with U and s from eigen decomp
template <typename T> 
void spectralDecomp<T>::initWithEigenDecomp( const T &U, const vec &s){

	this->U = U;
	this->s = s;
}

// Initialize with U and s from eigen decomp, and weights
template <typename T> 
void spectralDecomp<T>::initWithEigenDecomp( const T &U, const vec &s, const vec &weights){

	Rcpp::Rcout << "Reweighting not applied" << std::endl;

	// SVD after applying weights
	// mat Q, R;
	// qr(Q, R, scaleCols(U_, s_));

	this->U = U;
	this->s = s;
}




}




#endif