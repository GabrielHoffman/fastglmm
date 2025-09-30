
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
      s(other.s),
      Z(other.Z),
      Zw(other.Zw),
      V(other.V)
    {}

    /** Compute spectral decomp of weighted Z. Uses general method that applies to all (i.e. sparse, dense, discrete, continous) Z matrices.  
     * 
     * @param weights row weights
     * @param sort should singular values and vectors be sorted
    */ 
    void reweight( const vec &weights, const bool &sort = true){
      Rcpp::Rcout << "reweight spectralDecomp..." << std::endl;
      Zw = scaleEachCol(this->Z, sqrt(weights)); 
      eig_sym( this->s, V, mat(Zw.t() * Zw) );
      this->U = scaleEachRow(Zw * V, 1 / sqrt(this->s));

      if( sort ){
        // reorder by decreasing eigen-value
        uvec idx = sort_index(this->s, "descend");
        this->s = this->s(idx);
        this->U = this->U.cols(idx);
      }

      mat(Zw.t() * Zw).print("crossprod:");
      mat(U).print("U:");
      s.print("s:");
      mat(this->Z).print("Z:");
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
    T U;
    vec s;
    T Z, Zw;
    mat V;
};

/** Special case of spectralDecomp when Z is categorical
 */ 
template <typename T> 
class spectralDecompCategorical:
  public spectralDecomp<T> {
  public:

  spectralDecompCategorical(const T &U, const vec &s) :
    spectralDecomp<T>(U, s) {

    Rcpp::Rcout << "spectralDecompCategorical..."  << std::endl;
    mat(U).print("U:");
    mat(this->Z).print("Z:");
  }

  spectralDecompCategorical(const T &Z){
    this->Z = Z;
    vec ones(Z.n_rows, fill::ones);
    reweight(ones );
   }

  void reweight( const vec &weights, const bool &sort = true){
    Rcpp::Rcout << "reweight spectralDecompCategorical..." << std::endl;
    // SVD of new weighted Z
    this->s = (weights.t() * this->Z).t();
    this->U = scaleRowsCols(this->Z, sqrt(weights), 1 / sqrt(this->s));

    if( sort ){
      // reorder by decreasing eigen-value
      uvec idx = sort_index(this->s, "descend");
      this->s = this->s(idx);
      this->U = this->U.cols(idx);
    }
  }
};




}




#endif