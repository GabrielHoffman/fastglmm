/***************************************************************
 * @file    fastlmm_fit.h
 * @author  Gabriel Hoffman
 * @email   gabriel.hoffman@mssm.edu
 * @brief   Fit linear mixed model
 * Copyright (C) 2024 Gabriel Hoffman
 **************************************************************/

#ifndef _FASTLMM_FIT_H_
#define _FASTLMM_FIT_H_

// if -D USE_R, use RcppArmadillo library
#ifdef USE_R
// [[Rcpp::depends(RcppParallel)]]  
#include <RcppArmadillo.h>
#else
#include <armadillo>
#endif

using namespace arma;

#include "local_min.h"
#include "misc.h"
#include "ModelFit.h"
#include "spectralDecomp.h"
#include "satterthwaite.h"

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
            const spectralDecomp<T3> &dcmp,
            const vec &weights_,
            const ModelDetail md = LOW,
            const double &lambda = 0,
            const bool REML = false);

    // constructor, precompute Yu, Xu
    fastlmm(const T1 &Y_, 
            const T2 &X_, 
            const spectralDecomp<T3> &dcmp,
            const vec &weights_,
            const vec &Yu_, 
            const mat &Xu_,
            const ModelDetail md = LOW,
            const double &lambda = 0,
            const bool REML = false);

    // constructor, precompute Yu, Xu, Gamma_XX, Gamma_XY
    fastlmm(const T1 &Y_, 
            const T2 &X_,
            const spectralDecomp<T3> &dcmp,
            const vec &weights_,
            const vec &Yu_, 
            const mat &Xu_,
            const mat &Gamma_XX_, 
            const mat &Gamma_XY_,
            const ModelDetail md = LOW,
            const double &lambda = 0,
            const bool REML = false);

    // constructor without response
    fastlmm(const T2 &X_, 
            const spectralDecomp<T3> &dcmp,
            const ModelDetail md = LOW,
            const double &lambda = 0,
            const bool REML = false);

    void update_response(const T1 &Y_, const vec &weights_);
    void update_response(const T1 &Y_,
                         const vec &weights_, 
                         const mat &Yu_);

    // extract results  
    ModelFitLMM get_result(const bool &returnUS = false);

    // Accessors
    const double get_logLik() const { return this->logLik; }
    const vec get_beta() const { return this->beta; }
    const double get_sigSq_g() const { return sigSq_g;}
    const double get_sigSq_e() const { 
      return this->delta_hat * this->sigSq_g;
    }
    const int get_iter() const { return this->iter;}
    const double get_delta() const { return this->delta_hat;}
    
    const mat get_vcov() const {

      // initialize result matrix to NaN
      mat QXX_inv(this->QXX.n_rows, this->QXX.n_cols, arma::fill::nan);
      bool success = true;

      // if QXX is finite, compute inverse
      if( this->QXX.is_finite() ) {
        success = inv_sympd(QXX_inv, this->QXX);
      }

      // if inverse fails, set result to NaN
      if( !success ){
        QXX_inv.fill(arma::datum::nan);
      }

      // Scale NaN matrix
      return this->sigSq_g * QXX_inv ;
    }
    // condition number of vcov matrix
    const double get_vcov_kappa() const {
      return 1 / rcond(this->QXX) ;
    }
    const mat get_beta_se() const {
      return sqrt(diagvec(get_vcov()));
    }

    // if model fails, set beta to nan
    void set_model_failure(){ 
      beta.fill(datum::nan);
    }

    /** Residual degrees of freedom
     */ 
    const double get_rdf() const;

    const vec hatvalues() const; // diag of hat matrix
    const vec residuals() const; 
    const vec fitted() const;

    // Best linear unbiased predictor of random effect
    // same as ranef() in R
    const vec blup() const;

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
    T3 U, V;
    vec s, weights;
    T1 Yu;
    T2 Xu;
    mat Gamma_XX, Gamma_XY;
    vec inv_s_delta;
    mat inv_s_delta_Xu;
    mat QXX, QXY;
    vec beta;
    vec r, ru;
    ModelDetail md;
    double lambda;
    bool REML;
    double logLik, sigSq_g, delta_hat;
    int iter = 0;
    int n_active; // sample size with non-zero weight
    spectralDecomp<T3> dcmp;
};


// constructor, minimal
template <typename T1, typename T2, typename T3> 
fastlmm<T1, T2, T3>::fastlmm(
  const T1 &Y_, 
  const T2 &X_, 
  const spectralDecomp<T3> &dcmp,
  const vec &weights_,
  const ModelDetail md,
  const double &lambda,
  const bool REML):   
  wsqrt(sqrt(weights_)),
  Y(Y_ % wsqrt),
  X(scaleEachCol(X_, wsqrt)),
  weights(weights_),
  md(md), 
  lambda(lambda),
  REML(REML),
  dcmp(dcmp) {

  this->dcmp.reweight(weights);
  U = this->dcmp.get_U();
  s = this->dcmp.get_s();
  V = this->dcmp.get_V();

  n_active = accu(weights != 0.0);
  Yu = U.t() * Y;
  Xu = U.t() * X;
  Gamma_XX = X.t() * X - Xu.t() * Xu;
  Gamma_XY = X.t() * Y - Xu.t() * Yu; 
  inv_s_delta_Xu = mat( Xu.n_rows, Xu.n_cols);
} 




// constructor, precompute Yu, Xu
template <typename T1, typename T2, typename T3> 
fastlmm<T1, T2, T3>::fastlmm(
  const T1 &Y_, 
  const T2 &X_, 
  const spectralDecomp<T3> &dcmp,
  const vec &weights_,
  const vec &Yu_, 
  const mat &Xu_,
  const ModelDetail md,
  const double &lambda,
  const bool REML):   
  wsqrt(sqrt(weights_)),
  Y(Y_ % wsqrt),
  X(scaleEachCol(X_, wsqrt)),
  dcmp(dcmp),
  weights(weights_),
  md(md), 
  lambda(lambda),
  REML(REML) {

  this->dcmp.reweight(weights);
  U = this->dcmp.get_U();
  s = this->dcmp.get_s();
  V = this->dcmp.get_V();

  n_active = accu(weights != 0.0);
  Yu = Yu_;
  Xu = Xu_;
  Gamma_XX = X.t() * X - Xu.t() * Xu;
  Gamma_XY = X.t() * Y - Xu.t() * Yu; 
  inv_s_delta_Xu = mat( Xu.n_rows, Xu.n_cols);
} 

// constructor, precompute Yu, Xu, Gamma_XX, Gamma_XY
template <typename T1, typename T2, typename T3> 
fastlmm<T1, T2, T3>::fastlmm(
  const T1 &Y_, 
  const T2 &X_,
  const spectralDecomp<T3> &dcmp,
  const vec &weights_,
  const vec &Yu_, 
  const mat &Xu_,
  const mat &Gamma_XX_, 
  const mat &Gamma_XY_,
  const ModelDetail md,
  const double &lambda,
  const bool REML):   
  wsqrt(sqrt(weights_)),
  Y(Y_ % wsqrt),
  X(scaleEachCol(X_, wsqrt)),
  dcmp(dcmp),
  weights(weights_),
  md(md), 
  lambda(lambda),
  REML(REML) {

  this->dcmp.reweight(weights);
  U = this->dcmp.get_U();
  s = this->dcmp.get_s();
  V = this->dcmp.get_V();

  n_active = accu(weights != 0.0);
  Yu = Yu_;
  Xu = Xu_;
  Gamma_XX = Gamma_XX_;
  Gamma_XY = Gamma_XY_;
  inv_s_delta_Xu = mat( Xu.n_rows, Xu.n_cols);
} 


template <typename T1, typename T2, typename T3> 
  fastlmm<T1, T2, T3>::fastlmm( 
  const T2 &X_, 
  const spectralDecomp<T3> &dcmp,
  const ModelDetail md,
  const double &lambda,
  const bool REML): 
  X(X_),
  dcmp(dcmp),
  md(md), 
  lambda(lambda),
  REML(REML) { 

  this->dcmp.reweight(weights);
  U = this->dcmp.get_U();
  s = this->dcmp.get_s();
  V = this->dcmp.get_V();

  Xu = U.t() * X;
  Gamma_XX = X.t() * X - Xu.t() * Xu;
  inv_s_delta_Xu = mat( Xu.n_rows, Xu.n_cols);
}




template <typename T1, typename T2, typename T3> 
const double fastlmm<T1, T2, T3>::get_rdf() const {

  int n = n_active;
  int k = s.n_elem;

  // X is already scaled
  // sum(h1)
  // h1.sum <- with(object, delta*sum(1/(s+delta))) + (n-k)
  double h1_sum = delta_hat*(sum(1/(s+delta_hat)) + (n-k) / delta_hat);

  // sum(h2)
  // A <- with(object, X / delta - U %*% ((s/(delta*s + delta^2)) * crossprod(U, X)))
  // D <- solve(crossprod(A, X))
  // h2.sum <- object$delta * sum(A * (A %*% D))
  vec w = (s/(delta_hat*s + pow(delta_hat,2)));
  mat A = mat(X / delta_hat - U * scaleEachCol( Xu, w));
  double h2_sum = delta_hat * trace(solve(A.t() * X, A.t() * A));

  return h1_sum - h2_sum;
}

template <typename T1, typename T2, typename T3> 
const vec fastlmm<T1, T2, T3>::hatvalues() const {

  // Usq <- model$U^2
  T3 Usq = square(U);

  // h1 <- model$delta*with(model, Usq %*% (1/(s+delta))) + (1 - rowSums(Usq))
  vec h1 = delta_hat * Usq * (1/(s+delta_hat)) + (1 - sum(Usq, 1));

  // A <- with(model, X / delta - U %*% ((s/(delta*s + delta^2)) * crossprod(U, X)))
  // D <- solve(crossprod(A, X))
  // h2 <- model$delta * rowSums(A * (A %*% D))
  vec w = (s/(delta_hat*s + pow(delta_hat,2)));
  mat A = mat(X / delta_hat - U * scaleEachCol( Xu, w));
  mat D_A_t = solve(A.t() * X, A.t());
  vec h2 = delta_hat * sum(A % D_A_t.t(), 1);

  // hatvalues
  return 1 - h1 + h2;
}



template <typename T1, typename T2, typename T3> 
const vec fastlmm<T1, T2, T3>::residuals() const {

  return (Y / sqrt(weights)) - fitted();
}

// return predict(fit)
template <typename T1, typename T2, typename T3> 
const vec fastlmm<T1, T2, T3>::fitted() const {

  // ** need to scale X because it was transformed at the start
  // a <- object$U %*% (sqrt(object$s) * ranef.fastlmm(object))
  // a / sqrt(object$weights) + object$design %*% coef(object)
  return ((U * (sqrt(s) % blup())) + X * beta) / sqrt(weights);
}


template <typename T1, typename T2, typename T3> 
const vec fastlmm<T1, T2, T3>::blup() const {

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
  
  double n = n_active;
  double rank = Xu.n_rows;

  inv_s_delta = 1 / (s+delta);

  // inv_s_delta_Xu   <- inv_s_delta * Xu
  inv_s_delta_Xu = scaleEachCol(Xu, inv_s_delta);

  // QXX = crossprod(Xu, inv_s_delta_Xu) + Gamma_XX / delta
  QXX = Xu.t() * inv_s_delta_Xu + Gamma_XX / delta;

  // Ridge penalty
  // QXX.diag() += lambda;
  // but not on intercept
  for( uword i=1; i<QXX.n_rows; ++i){
    QXX(i,i) += lambda;
  }

  // QXY = crossprod(Xu, inv_s_delta_Yu) + Gamma_XY / delta
  QXY = Xu.t() * (inv_s_delta % Yu) + Gamma_XY / delta;

  // beta <<- solve( QXX, QXY)
  // beta = solve(QXX, QXY, solve_opts::likely_sympd);
  int status = solve(beta, QXX, QXY, 
    solve_opts::likely_sympd + arma::solve_opts::no_approx);

  // if model failed
  if( ! status ){
    // set beta to NAN
    beta.set_size(QXX.n_rows);
    beta.fill(datum::nan);
  }

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
  this->logLik = -1*local_min(leftIn, rightIn, tol, &F, res, iter);

  // augment with value this is constant for varying delta's
  // weights with zero value, give Inf log values
  // so use omit_nonfinite
  this->logLik += sum(omit_nonfinite(log(weights)))/2.0;

  this->delta_hat = exp(res);
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

  n_active = accu(weights_ != 0.0);
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
                      1.0,
                      get_beta());

  // set additional values based on ModelDetail md
  mat V = this->get_vcov();

  switch( md ){
    case MAX:       
    case MOST:
      res.hatvalues = this->hatvalues(); 
    case HIGH: 
      res.residuals = this->residuals();
    case MEDIUM: 
      res.vcov = V;
    case LOW: 
      res.se = sqrt(diagvec(V)); 
      res.rdf = this->get_rdf();
    case LEAST: 
      break;
  }

  // Precompute values for Satterthwaite DDF to be
  // used later with V and L specified
  Satterthwaite ddf_sat(res.y.n_elem, res.sigSq_g, res.sigSq_e, s, Xu, Gamma_XX, inv_s_delta, inv_s_delta_Xu);

  // save precomputed values
  res.hessian_vc = ddf_sat.get_hessian();
  res.A_sat = ddf_sat.get_A();
  res.B_sat = ddf_sat.get_B();

  // if returnUS
  // return U and s 
  if( returnUS ){
    switch( dcmp.get_type() ){
      case GENERAL:
        res.setUS(U, s, dcmp.get_V());
        break;
      case CATEGORICAL:
        // V is identity
        res.setUS(U, s, eye<sp_mat>(U.n_cols, U.n_cols));
        break;
    };
  }

  return res;
}



} // end namespace


#endif
