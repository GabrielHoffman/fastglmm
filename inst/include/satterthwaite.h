/***************************************************************
 * @file    satterthwaite.h
 * @author  Gabriel Hoffman
 * @email   gabriel.hoffman@mssm.edu
 * @brief   Precompute values for Satterthwaite approx to ddf
 * Copyright (C) 2024 Gabriel Hoffman
 **************************************************************/

#ifndef _SATTERTHWAITE_H_
#define _SATTERTHWAITE_H_

// if -D USE_R, use RcppArmadillo library
#ifdef USE_R
// [[Rcpp::depends(RcppParallel)]]  
#include <RcppArmadillo.h>
#else
#include <armadillo>
#endif

using namespace arma;

namespace fastglmmLib {

// denominator degrees of freedom using Satterthwaite method
// For math, see https://chatgpt.com/share/698f66b1-88b0-800b-9001-e28f444f6015  
template <typename T> 
class Satterthwaite {
  public:
    // constructor
    Satterthwaite( 
      const int &n, 
      const double &sigSq_g,
      const double &sigSq_e,
      const vec &s, 
      const T &Xu, 
      const mat &Gamma_XX, 
      const vec &inv_s_delta,
      const mat &inv_s_delta_Xu);

    /** Evalute denominator degrees of freedom
     */
    vec get_ddf(      
      const mat &V,
      const mat &L) const;

    /** Get precomputed values to compute DDF from V and L later
     */ 
    const mat get_hessian() const {return H;}
    const mat get_A() const {return A;}
    const mat get_B() const {return B;}

  private:
    /* Hessian of sigma^2_g, sigma^2_e
    */
    void set_hessian(
      const vec &s);

    /* Precompute value for gradient
    */
    void prep_gradient(
      const vec &s, 
      const T &Xu, 
      const vec &inv_s_delta,
      const mat &inv_s_delta_Xu,
      const mat &Gamma_XX);

    /** Evalute gradient
     */
    vec get_gradient( 
      const mat &L) const;

    int n;
    double sigSq_g, sigSq_e, delta;
    mat A, B, H;
};


// constructor
template <typename T> 
Satterthwaite<T>::Satterthwaite( 
  const int &n, 
  const double &sigSq_g,
  const double &sigSq_e,    
  const vec &s, 
  const T &Xu, 
  const mat &Gamma_XX, 
  const vec &inv_s_delta,
  const mat &inv_s_delta_Xu):
  n(n),
  sigSq_g(sigSq_g),
  sigSq_e(sigSq_e),
  delta(sigSq_e/sigSq_g) {

  set_hessian(s);
  prep_gradient(s, Xu, inv_s_delta, inv_s_delta_Xu, Gamma_XX);
}


// denominator degrees of freedom using Satterthwaite method 
template <typename T> 
vec Satterthwaite<T>::get_ddf( const mat &V, const mat &L) const{

  mat H_inv = inv(H);
  mat g = get_gradient(L);

  vec df(L.n_rows);
  double var_Lbeta, v_numerator, v_denom;

  // for each contrast (i.e. row)
  for(int i=0; i<L.n_rows; i++){
    // var_Lbeta <- crossprod(L[i,], vcov(fit)) %*% L[i,]
    var_Lbeta = as_scalar(L.row(i).t() * V * L.row(i));
    v_numerator = 2 * pow(var_Lbeta, 2);

    // v_denom <- crossprod(g[[i]], A) %*% g[[i]]
    v_denom = as_scalar(g.row(i).t() * H_inv * g.row(i));

    df(i) = v_numerator / v_denom;
  }

  return df;
}

template <typename T> 
vec Satterthwaite<T>::get_gradient(const mat &L) const {

  vec g(L.n_rows, 2, fill::zeros);
  mat invAL;
  double c;

  // for each contrast (i.e. row)
  for(int i=0; i<L.n_rows; i++){
    // invAL <- solve(A, L[i,])
    invAL = solve(A, L.row(i));

    // C <- invAL %*% B %*% invAL
    c = as_scalar(invAL * B * invAL);

    // crossprod(L[i,], invAL) - fit$delta*C
    g(i,0) = as_scalar(L.row(i).t() * invAL - delta*c);
    g(i,1) = c;
  }

  return g;
}

template <typename T> 
void Satterthwaite<T>::set_hessian(const vec &s){

  // n <- length(fit$y)
  // r <- length(s)
  double r = s.n_elem;

  // A <- sum(1/(s + delta)) + (n-r)/delta
  double a = accu(1/(s + delta)) + (n-r)/delta;

  // B <- sum(1/(s + delta)^2) + (n-r)/delta^2
  double b = accu(1/pow(s + delta,2)) + (n-r)/pow(delta, 2);

  H = mat(2,2, fill::zeros);

  // H[1,1] <- (n - 2*delta*a + delta^2*b) / (2*sigSq_g^2)
  H(0,0) = (n - 2*delta*a + pow(delta,2)*b) / (2*pow(sigSq_g,2));

  // H[1,2] <- H[2,1] <- (a - delta*b) / (2*sigSq_g^2)
  H(0,1) = (a - delta*b) / (2*pow(sigSq_g,2));
  H(1,0) = H(0,1);

  // H[2,2] <- b / (2*sigSq_g^2)
  H(1,1) = b / (2*pow(sigSq_g,2));
}



// gradient of variance for Satterthwaite method

template <typename T> 
void Satterthwaite<T>::prep_gradient(
  const vec &s, 
  const T &Xu, 
  const vec &inv_s_delta,
  const mat &inv_s_delta_Xu,
  const mat &Gamma_XX) {

  // Xu <- crossprod(fit$U, fit$design)
  // Gamma_XX <- crossprod(X) - crossprod(Xu)
  // inv_s_delta <- 1 / (fit$s + fit$delta)
  // inv_s_delta_Xu <- inv_s_delta * Xu

  // A <- crossprod(Xu, inv_s_delta_Xu) + Gamma_XX / fit$delta
  // inv_s_delta_Xu = scaleEachCol(Xu, inv_s_delta);
  A = Xu.t() * inv_s_delta_Xu + Gamma_XX / delta;

  // inv_s_delta_Xu <- inv_s_delta^2 * Xu
  // B <- crossprod(Xu, inv_s_delta_Xu) + Gamma_XX / fit$delta^2
  B = Xu.t() * scaleEachCol(Xu, pow(inv_s_delta,2)) + Gamma_XX / pow(delta, 2);
}





} // end namespace


#endif
