

#ifndef _GLMM_LAPLACE_
#define _GLMM_LAPLACE_

#ifdef USE_R
// [[Rcpp::depends(RcppParallel, roptim)]]  
#include <RcppArmadillo.h>
#include <roptim.h>
#else
#include <armadillo>
#endif

#include "distributions.h"

using namespace arma;
using namespace std;
using namespace fastglmmLib;
using namespace roptim;

namespace fastglmmLib {

// Order of template variables
// T1 y
// T2 X
// T3 U
template <typename T1, typename T2, typename T3> 
class GLMMLaplace : public Functor {

  public:

  // Constructor
  GLMMLaplace( 
    const T1 &y, 
    const T2 &X, 
    const spectralDecomp<T3> &dcmp,
    const vec &weights,
    const vec &offset,
    const double &lambda,
    const string &family, 
    const ModelDetail md = LOW,
    GLMWork *work = NULL):
  y(y), 
  X(X),
  offset(offset), 
  dcmp(dcmp), 
  family(family), 
  work(work),
  lambda(lambda)
  {
    Z = (dcmp.get_U() * scaleEachCol(dcmp.get_V(), sqrt(dcmp.get_s())));

    // initialize to zero
    u = vec(Z.n_cols);
    u.fill(0.0); 
  }

  const vec get_u(){
    return u;
  }
  const T3 get_Z(){
    return Z;
  }

  void set_theta(const double &value){
    theta_fixed = value;
  }

  // negative log-likelihood
  double operator()(const arma::vec &x) override {

    vec beta = x.subvec(span(0, X.n_cols-1));
    double sig_b = exp(x[X.n_cols]);

    if( theta_fixed != -1){
      theta = theta_fixed;
    }else{
      // value of theta is estimated
      theta = exp(x[X.n_cols+1]);
    }

    auto [u_local, mu, HlogDet] = mode_u(beta, sig_b, theta);
    this->u = u_local;

    // ll_cond <- sum(
      // dnbinom(y, mu = res$mu, size = theta, log = TRUE)) 
    vec prob = theta/(theta+mu);
    double ll_cond = 0;
    for(int i=0; i<y.n_elem; i++){
      ll_cond += dnbinom_boost(y[i], theta, prob[i], true);
    }

    // ll_prior <- -0.5 * crossprod(res$u)[1] #sum(fit$u^2)
    double ll_prior = -0.5*dot(u_local,u_local);

    // ridge prior
    ll_prior += - 0.5*lambda*dot(beta, beta);

    // ll_det <- -0.5 * res$H.logDet 
    double ll_det = -0.5*HlogDet;

    // return the negative log-likelihood
    return -1*(ll_cond + ll_prior + ll_det);
  }

  private:
  vec y;
  mat X;
  vec offset;
  spectralDecomp<T3> dcmp;
  string family;
  double lambda;
  GLMWork *work;
  T3 Z;
  vec u;
  vec rhs, u_new, Xb_offset;
  sp_mat H;
  double theta_fixed = -1;
  double theta;

  // Compute mode of random effect u 
  tuple<vec, vec, double> mode_u(
    const vec &beta,
    const double &sig_b, 
    const double &theta,
    const double tol = 1e-4,
    const int maxit = 500)
  {

    // set to zero 
    // u.fill(0.0); 

    T3 ZL = sig_b*Z;
    Xb_offset = X * beta + offset; 

    for( int iter=0; iter<maxit; iter++){

      // eta <- drop(offset + XBeta + ZL %*% u)
      // eta <- pmax(pmin(eta, 12), -12) # clamp
      // mu  <- exp(eta)
      work->eta = Xb_offset + ZL * u;
      work->eta.clamp(-25, 25); // avoid underflow of mu
      work->mu = exp(work->eta);

      // work->w = work->mu / (1.0 + work->mu/theta);
      work->w = theta * work->mu / (theta + work->mu);

      work->z = work->eta + (y - work->mu) / work->mu;

      // H <- crossprod(ZL, ZL * w) + Diagonal(q, 1)
      H = ZL.t() * scaleEachCol(ZL, work->w);
      H.diag() += 1;

      // rhs <- crossprod(ZL, w * (z - XBeta - offset))
      rhs = ZL.t() * (work->w % (work->z - Xb_offset));

      // u_new <- solve(H, rhs) 
      // u_new = spsolve(H, rhs);
      // Since Z is categorical H is diagonal
      u_new = rhs / H.diag(); 

      if (max(abs(u_new - u)) < tol) break;

      u = u_new;
    }

    work->eta = offset + X * beta + ZL * u;
    work->eta.clamp(-25, 25); // avoid underflow of mu
    work->mu = exp(work->eta);
    double HlogDet = log_det_sympd(mat(H));

    return {u, work->mu, HlogDet};
  }

};


} // end namespace
#endif
