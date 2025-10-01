
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

typedef enum {
  GENERAL,
  CATEGORICAL
} ZTYPE;

/**
 * Class storing matrix Z, its _squared_ singular values (s), left singular vectors (U), and right singular vectors (V).  Z can be mat or sp_mat.  Default type is GENERAL for arbitrary matrix Z.  If Z is a categorical design matrix, use a faster method with type CATEGORICAL. 

 reweight() computes the SVD of the weighted Z matrix
 */ 
template <typename T> 
class spectralDecomp {

  public:
    spectralDecomp() {};

    /** Constructor 
     * @param U left singular vectors
     * @param s squared singular values
     */
    spectralDecomp(
      const T &U, 
      const vec &s, 
      const ZTYPE &type = GENERAL) : 
      U(U), 
      s(s), 
      Z(scaleEachRow(U, sqrt(s))),
      type(type) 
      {}

    /** Constructor 
     * @param Z Random effects design matrix
     */
    spectralDecomp( 
      const T &Z, 
      const ZTYPE &type = GENERAL ) : 
      Z(Z),      
      type(type){

      vec one(Z.n_rows, fill::ones);
      reweight(one);
    }

    // Copy constructor
    spectralDecomp(const spectralDecomp& other):
      U(other.U),
      s(other.s),
      Z(other.Z),
      Zw(other.Zw),
      V(other.V),
      type(other.type){
    }

    /** Compute spectral decomp of weighted Z. 
     * 
     * @param weights row weights
     * @param sort should singular values and vectors be sorted
    */ 
    void reweight( const vec &weights, const bool &sort = false){

      switch( type ){
        case GENERAL:
          Zw = scaleEachCol(Z, sqrt(weights)); 
          eig_sym( s, V, mat(Zw.t() * Zw) );
          U = scaleEachRow(Zw * V, 1 / sqrt(s));
          break;
        case CATEGORICAL:
          s = (weights.t() * Z).t();
          U = scaleRowsCols(Z, sqrt(weights), 1 / sqrt(s));
          break;
      }

      // mat Z_recon = scaleRowsCols(U, 1 / sqrt(weights), sqrt(s)) * V.t();
      // Z_recon.print("Z_recon:");

      if( sort ){
        // reorder by decreasing eigen-value
        uvec idx = sort_index(s, "descend");
        s = s(idx);
        U = U.cols(idx);
      }
    }

    /** Accessor, returns left singular vectors */
    T get_U() const {
      return U;
    }

    /** Accessor, returns squared singular values */
    vec get_s() const {
      return s;
    }

    /** Accessor, returns right singular values */
    mat get_V() const {

      mat V_ret;

      switch( type ){
        case GENERAL:
          V_ret = V;
          break;
        case CATEGORICAL:      
          V_ret = eye<mat>(U.n_cols, U.n_cols);
          break;
      }

      return V_ret;
    }

    ZTYPE get_type() const {
      return type;
    }

  protected:
    T U;
    vec s;
    T Z, Zw;
    mat V;
    ZTYPE type;
};




}




#endif