
#ifndef SPECTRAL_DECOMP_H_
#define SPECTRAL_DECOMP_H_

#include <armadillo>
using namespace arma;

#include "misc.h"

namespace fastlmmLib {


template <typename T> 
class spectralDecomp {

  public:
  	spectralDecomp(){};

  	// Initialize with indicator matrix Z;
    void initWithIndicator( const T &Z_, const vec &weights_);

    // Initialize with U and s from eigen decomp
    void initWithEigenDecomp( const T &U_, const vec &s_);

    // Initialize with U and s from eigen decomp, and weights
    void initWithEigenDecomp( const T &U_, const vec &s_, const vec &weights_);

    T get_vectors(){return U;}
    vec get_values(){return s;}

  private:
   T U;
   vec s;
};

// Initialize with indicator matrix
template <typename T> 
void spectralDecomp<T>::initWithIndicator( const T &Z_, const vec &weights_){

	s = (weights_.t() * Z_).t();
	U = scaleRowsCols(Z_, sqrt(weights_), 1 / sqrt(s));
}


// Initialize with U and s from eigen decomp
template <typename T> 
void spectralDecomp<T>::initWithEigenDecomp( const T &U_, const vec &s_){

	this->U = U_;
	this->s = s_;
}

// Initialize with U and s from eigen decomp, and weights
template <typename T> 
void spectralDecomp<T>::initWithEigenDecomp( const T &U_, const vec &s_, const vec &weights_){

	// SVD after applying weights
	// mat Q, R;
	// qr(Q, R, scaleCols(U_, s_));

	// this->U = U_;
	// this->s = s_;
}




}




#endif