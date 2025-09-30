
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

/**
 * Class storing matrix Z, its squared singular values (s) and the left singular vectors (U). Z can be mat or sp_mat.  spectralDecomp.reweight(w) uses a general method that applies to all (i.e. sparse, dense, discrete, continous) Z matrices.  If Z is sparse, a faster method can be applied using the spectralDecompDiscrete class.  
 */ 
template <typename T> 
class spectralDecomp {

  public:
    spectralDecomp() {};

    /** Constructor 
     * @param U left singular vectors
     * @param s squared singular values
     */
    spectralDecomp(const T &U, const vec &s) : 
      U(U), s(s), Z(scaleEachRow(U, sqrt(s))) 
      {}

    /** Constructor 
     * @param Z Random effects design matrix
     */
    spectralDecomp( const T &Z ) : 
      Z(Z) {

      vec one(Z.n_rows, fill::ones);
      reweight(one);
    }

    // Copy constructor
    spectralDecomp(const spectralDecomp& other): 
      U(other.U),
      Z(other.Z),
      Zw(other.Zw),
      V(other.V),
      s(other.s)
    {}

    /** Compute spectral decomp of weighted Z. Uses general method that applies to all (i.e. sparse, dense, discrete, continous) Z matrices.  
     * 
     * @param weights row weights
    */ 
    void reweight( const vec &weights){
      Rcpp::Rcout << "reweight spectralDecomp..." << std::endl;
      Zw = scaleEachCol(this->Z, sqrt(weights)); 
      eig_sym( this->s, V, mat(Zw.t() * Zw) );
      this->U = scaleEachRow(Zw * V, 1 / sqrt(this->s));
    }

    /** Accessor, returns left singular vectors */
    T get_vectors() const {
      return U;
    }

    /** Accessor, returns squared singular values */
    vec get_values() const {
      return s;
    }

  protected:
    T U, Z, Zw;
    mat V;
    vec s;
};

/** Special case of spectralDecomp when Z is categorical
 */ 
template <typename T> 
class spectralDecompCategorical:
  public spectralDecomp<T> {
  public:

  spectralDecompCategorical(const T &U, const vec &s) :
    spectralDecomp<T>(U, s) {}

  spectralDecompCategorical(const T &Z){
    this->Z = Z;
    vec ones(Z.n_rows, fill::ones);
    reweight(ones );
   }

  void reweight( const vec &weights){
    Rcpp::Rcout << "reweight spectralDecompCategorical..." << std::endl;
    // SVD of new weighted Z
    this->s = (weights.t() * this->Z).t();
    this->U = scaleRowsCols(this->Z, sqrt(weights), 1 / sqrt(this->s));
  }
};




}




#endif