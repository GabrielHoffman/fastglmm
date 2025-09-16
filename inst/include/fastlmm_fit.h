/***************************************************************
 * @file    fastlmm_fit.h
 * @author  Gabriel Hoffman
 * @email   gabriel.hoffman@mssm.edu
 * @brief   Fit linear mixed model
 * Copyright (C) 2024 Gabriel Hoffman
 **************************************************************/

#ifndef _FASTLMM_FIT_H_
#define _FASTLMM_FIT_H_

// if -D ARMA, use plain armadillo library
#ifdef ARMA
#include <armadillo>
#else
#include <RcppArmadillo.h>
#endif

using namespace arma;

#include "local_min.h"
#include "misc.h"
#include "ModelFit.h"

namespace fastglmmLib {

// Order of template variables
// T1 Y
// T2 X
// T3 U
template <typename T1, typename T2, typename T3> 
class fastlmm {       
  public:  

    // constructor, minimal
    fastlmm(){};

    fastlmm(const T1 &Y_, 
            const T2 &X_, 
            const T3 &U_, 
            const vec &s_,
            const vec &weights_,
            const ModelDetail md = LOW,
            const bool REML = false);

    // constructor, precompute Yu, Xu
    fastlmm(const T1 &Y_, 
            const T2 &X_, 
            const T3 &U_, 
            const vec &s_,
            const vec &weights_,
            const vec &Yu_, 
            const mat &Xu_,
            const ModelDetail md = LOW,
            const bool REML = false);

    // constructor, precompute Yu, Xu, Gamma_XX, Gamma_XY
    fastlmm(const T1 &Y_, 
            const T2 &X_,
            const T3 &U_, 
            const vec &s_,
            const vec &weights_,
            const vec &Yu_, 
            const mat &Xu_,
            const mat &Gamma_XX_, 
            const mat &Gamma_XY_,
            const ModelDetail md = LOW,
            const bool REML = false);

    // constructor without response
    fastlmm(const T2 &X_, 
            const T3 &U_, 
            const vec &s_,
            const ModelDetail md = LOW,
            const bool REML = false);

    void update_response(const T1 &Y_, const vec &weights_);
    void update_response(const T1 &Y_,
                         const vec &weights_, 
                         const mat &Yu_);

    // extract results  
    ModelFitLMM get_result(const bool &returnUS = false);

    // Accessors
    const double get_logLik(){ return this->logLik; }
    const vec get_beta(){ return this->beta; }
    const double get_sigSq_g(){ return sigSq_g;}
    const double get_sigSq_e(){ 
      return this->delta_hat * this->sigSq_g;
    }
    const int get_iter(){ return this->iter;}
    const double get_delta(){ return this->delta_hat;}
    const mat get_vcov(){
      return inv_sympd(this->QXX) * this->sigSq_g;
    }
    const mat get_beta_se(){
      return sqrt(diagvec(get_vcov()));
    }

    // how to combine X and K?
    // to create diagonals of hat matrix?
    const double get_edf(); // defined before
    const double get_rdf(){
      return r.n_rows - X.n_cols;
    } // based on Hastie, et al
    const vec hatvalues(); // diag of hat matrix
    const vec residuals();
    const vec fitted();
    // df = sum(s[seq_len(rank)]/(s[seq_len(rank)]+delta))

    // Best linear unbiased predictor of random effect
    // same as ranef() in R
    const vec blup();

    // compute log likelihood
    double ll(const double &delta);

    void estimate_delta(  const double &left,
                          const double &right,
                          const double &tol);
    // evaluate logLik, beta, etc at delta value
    void eval_delta( const double &delta){
      this->logLik = ll( delta );
      this->delta_hat = delta;
    }

    // Score test
    double score_test(const vec &x_);

    // Update Y, keeping rest constant
    void update_Y( const T1 &Y_);

    // Update X, keeping rest constant
    void update_X( const vec &X_);

    vec get_weights(){ return weights;}

    vec get_ru(){ return ru;}
    vec get_r(){ return r;}
    vec get_y(){ return Y;}

  private:
    vec wsqrt;
    T1 Y;
    T2 X;
    T3 U;
    vec s, weights;
    T1 Yu;
    T2 Xu;
    mat Gamma_XX, Gamma_XY;
    vec inv_s_delta;
    mat inv_s_delta_Xu;
    mat QXX, QXY;
    mat beta;
    vec r, ru;
    ModelDetail md;
    bool REML;
    double logLik, sigSq_g, delta_hat;
    int iter = 0;
};


// constructor, minimal
template <typename T1, typename T2, typename T3> 
fastlmm<T1, T2, T3>::fastlmm(const T1 &Y_, 
        const T2 &X_, 
        const T3 &U_, 
        const vec &s_,
        const vec &weights_,
        const ModelDetail md,
        const bool REML):   
  wsqrt(sqrt(weights_)),
  Y(Y_ % wsqrt),
  X(scaleEachCol(X_, wsqrt)),
  U(U_),
  s(s_),
  weights(weights_),
  md(md), 
  REML(REML) {
  Yu = U_.t() * Y;
  Xu = U_.t() * X;
  Gamma_XX = X.t() * X - Xu.t() * Xu;
  Gamma_XY = X.t() * Y - Xu.t() * Yu; 
  inv_s_delta_Xu = mat( Xu.n_rows, Xu.n_cols);
} 




// constructor, precompute Yu, Xu
template <typename T1, typename T2, typename T3> 
fastlmm<T1, T2, T3>::fastlmm(const T1 &Y_, 
        const T2 &X_, 
        const T3 &U_, 
        const vec &s_,
        const vec &weights_,
        const vec &Yu_, 
        const mat &Xu_,
        const ModelDetail md,
        const bool REML):   
  wsqrt(sqrt(weights_)),
  Y(Y_ % wsqrt),
  X(scaleEachCol(X_, wsqrt)),
  U(U_),
  s(s_),
  weights(weights_),
  md(md), 
  REML(REML) {
  Yu = Yu_;
  Xu = Xu_;
  Gamma_XX = X.t() * X - Xu.t() * Xu;
  Gamma_XY = X.t() * Y - Xu.t() * Yu; 
  inv_s_delta_Xu = mat( Xu.n_rows, Xu.n_cols);
} 

// constructor, precompute Yu, Xu, Gamma_XX, Gamma_XY
template <typename T1, typename T2, typename T3> 
fastlmm<T1, T2, T3>::fastlmm(const T1 &Y_, 
                            const T2 &X_,
                            const T3 &U_, 
                            const vec &s_,
                            const vec &weights_,
                            const vec &Yu_, 
                            const mat &Xu_,
                            const mat &Gamma_XX_, 
                            const mat &Gamma_XY_,
                            const ModelDetail md,
                            const bool REML):   
  wsqrt(sqrt(weights_)),
  Y(Y_ % wsqrt),
  X(scaleEachCol(X_, wsqrt)),
  U(U_),
  s(s_),
  weights(weights_),
  md(md), 
  REML(REML) {
  Yu = Yu_;
  Xu = Xu_;
  Gamma_XX = Gamma_XX_;
  Gamma_XY = Gamma_XY_;
  inv_s_delta_Xu = mat( Xu.n_rows, Xu.n_cols);
} 


template <typename T1, typename T2, typename T3> 
fastlmm<T1, T2, T3>::fastlmm( const T2 &X_, 
                              const T3 &U_, 
                              const vec &s_,
                              const ModelDetail md,
                              const bool REML): md(md), REML(REML) {
  X = X_;
  U = U_;
  s = s_;
  Xu = U_.t() * X_;
  Gamma_XX = X.t() * X - Xu.t() * Xu;
  inv_s_delta_Xu = mat( Xu.n_rows, Xu.n_cols);
}

template <typename T1, typename T2, typename T3> 
const vec fastlmm<T1, T2, T3>::hatvalues(){
  vec h(Y.n_elem, fill::ones);
  return h;
}

template <typename T1, typename T2, typename T3> 
const vec fastlmm<T1, T2, T3>::residuals(){
  return (Y / sqrt(weights)) - fitted();
}


template <typename T1, typename T2, typename T3> 
const vec fastlmm<T1, T2, T3>::fitted(){

  // a <- object$U %*% (sqrt(object$s) * ranef.fastlmm(object))
  // a / sqrt(object$weights) + object$design %*% coef(object)
  // need to scale X because it was transformed at the start
  return ((U * (sqrt(s) % blup())) + X * beta) / sqrt(weights);
}


template <typename T1, typename T2, typename T3> 
const vec fastlmm<T1, T2, T3>::blup(){

  // Zw <- c(sqrt(fit$weights)) * fit$Z
  // A <- crossprod(fit$U, Zw)
  // A <- with(fit, crossprod(fit$U, c(sqrt(weights)) * U * sqrt(s)))
  // b <- fit$ru / (fit$s + fit$delta)
  // crossprod(A, b)

  // T3 A = U.t() * scaleRowsCols(U, sqrt(weights), sqrt(s));
  // vec b = ru / (s + delta_hat);
  // return A.t() * b;

  // # since U^T U is identity if the GRM is full rank
  // v <- with(object, sqrt(s)*ru / (s + delta))

  return (sqrt(s) % ru) / (s + delta_hat);
}




template <typename T1, typename T2, typename T3>  
double fastlmm<T1, T2, T3>::ll(const double &delta ) { 

  double n = X.n_rows;
  double rank = Xu.n_rows;

  inv_s_delta = 1 / (s+delta);

  // inv_s_delta_Xu   <- inv_s_delta * Xu
  inv_s_delta_Xu = scaleEachCol(Xu, inv_s_delta);

  // QXX = crossprod(Xu, inv_s_delta_Xu) + Gamma_XX / delta
  QXX = Xu.t() * inv_s_delta_Xu + Gamma_XX / delta;

  // QXY = crossprod(Xu, inv_s_delta_Yu) + Gamma_XY / delta
  QXY = Xu.t() * (inv_s_delta % Yu) + Gamma_XY / delta;

  // beta <<- solve( QXX, QXY)
  beta = solve(QXX, QXY, solve_opts::likely_sympd);

  // # Eval sig_g
  // ru <- Yu - Xu %*% beta
  ru = Yu - Xu * beta;

  // r <- Y - X %*% beta
  r = Y - X * beta;

  // Qrr <- crossprod(ru, inv_s_delta_ru) + (crossprod(r)[1] - crossprod(ru)[1])/ delta
  // sig_g <<- Qrr[1] / n
  double QRR = dot(ru, (inv_s_delta % ru)) + (dot(r,r) - dot(ru,ru)) / delta;
  sigSq_g = QRR / n;

  // use 2.0 to ensure double precision
  double logLik = -n/2.0 * log(2.0*M_PI*sigSq_g) - 1.0/2.0 * (sum( log(s + delta ) ) + (n-rank) * log(delta)) - n/2.0; 

  // this is fixed, so don't eval every time, 
  //     just after estimation
  // + sum(log(weights))/2.0;

  return logLik;
}
    


// function to be minimized
static inline double ll_alone_mat( double delta_log, void *arg){

  auto *fit = (fastlmm<mat,mat,mat> *) arg;

  // search is done in log space 
  //  to give faster convergence
  fit->eval_delta( exp(delta_log) );

  return -1.0*fit->get_logLik();
}

// sparse version
static inline double ll_alone_spmat( double delta_log, void *arg){

  auto *fit = (fastlmm<mat,mat,sp_mat> *) arg;

  // search is done in log space 
  //  to give faster convergence
  fit->eval_delta( exp(delta_log) );

  return -1.0*fit->get_logLik();
}


template <typename T1, typename T2, typename T3> 
void fastlmm<T1, T2, T3>::estimate_delta( const double &left, const double &right, const double &tol ){

  double leftIn = left; 
  double rightIn = right;
  iter = 0;
  
  // initialize function
  funcStruct F;  
  F.params = this;

  // Since F.function can't take templated function
  if( isSpMatrix( U ) ){
    F.function = & ll_alone_spmat;
  }else{
    F.function = & ll_alone_mat;
  }

  // get maximize log-likelihood
  // need to mutliply but -1 since it actually minimizes
  // evaluated at minimum value 
  double res;
  logLik = -1*local_min(leftIn, rightIn, tol, &F, res, iter);

  // augment with value this is constant for varying delta's
  logLik += sum(log(weights))/2.0;

  delta_hat = exp(res);
}


template <typename T1, typename T2, typename T3> 
void fastlmm<T1, T2, T3>::update_response(const T1 &Y_,
                                          const vec &weights_){

  update_response(Y, weights_, U.t() * Y_);
} 


template <typename T1, typename T2, typename T3> 
void fastlmm<T1, T2, T3>::update_response(const T1 &Y_,
                                          const vec &weights_,
                                          const mat &Yu_){

  // indicator_decomp
  // modiy this->U  and this->s internally
  // compute sqrt(weights) for 
  // Y <- Y * sqrt(weights)
  // X <- X * sqrt(weights)
  // vec sqrtW = sqrt(weights_);
  // update_weights( Y_, X_, U_, s_, weights_);
  // Need to save X, U, s unmodified so it
  // can be weighted later


  this->weights = weights_;
  this->Y = Y_;
  this->Yu = Yu_;  
  this->Gamma_XY = X.t() * Y - Xu.t() * Yu;
}

 
template <typename T1, typename T2, typename T3> 
ModelFitLMM fastlmm<T1, T2, T3>::get_result(
          const bool &returnUS){

  // initialize with standard entries
  ModelFitLMM res = ModelFitLMM( true, 
                      get_logLik(),
                      get_weights(),
                      get_ru(),
                      get_y(),
                      get_delta(),
                      get_sigSq_g(),
                      get_sigSq_e(),
                      get_iter(),
                      get_beta());

  // res.dispersion = get_sigSq_e();

  // set additional values based on ModelDetail md
  mat V = get_vcov();

  switch( md ){
    case MAX: 
    case MOST:
      res.hatvalues = hatvalues(); 
    case HIGH: 
      res.residuals = residuals();
    case MEDIUM: 
      res.vcov = V;
    case LOW: 
      res.se = sqrt(diagvec(V)); 
      res.rdf = get_rdf();
    case LEAST: 
      break;
  }

  // if returnUS
  // return U and s 
  if( returnUS ){
    res.setUS(U, s);
  }

  return res;
}



} // end namespace


#endif
