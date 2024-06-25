#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifndef FASTLMM_H_
#define FASTLMM_H_

#include <RcppGSL.h>

#include <gsl/gsl_min.h>
#include <gsl/gsl_errno.h>

using namespace Rcpp; 
using namespace arma;

class FASTLMM {       
  public:  
    // constructor, minimal
    FASTLMM(const arma::vec &Y_, 
            const arma::mat &X_, 
            const arma::mat &U_, 
            const arma::vec &s_){

      Y = Y_;
      X = X_;
      U = U_;
      s = s_;

      // use weights here
      Yu = U.t() * Y;
      Xu = U.t() * X;

      // cp_X_low = crossprod(X) - crossprod(Xu)
      cp_X_low = X.t() * X - Xu.t() * Xu; 

      // cp_X_low_Y_low = crossprod(X, Y) - crossprod(Xu, Yu)
      cp_X_low_Y_low = X.t() * Y - Xu.t() * Yu; 

      inv_s_delta_Xu = arma::mat( Xu.n_rows, Xu.n_cols);
    } 

    // constructor, precompute Yu, Xu
    FASTLMM(const arma::vec &Y_, 
            const arma::mat &X_, 
            const arma::vec &Yu_, 
            const arma::mat &Xu_,
            const arma::mat &U_, 
            const arma::vec &s_){

      Y = Y_;
      X = X_;
      U = U_;
      s = s_;
      Yu = Yu_;
      Xu = Xu_;

      // use weights here

      // cp_X_low = crossprod(X) - crossprod(Xu)
      cp_X_low = X.t() * X - Xu.t() * Xu; 

      // cp_X_low_Y_low = crossprod(X, Y) - crossprod(Xu, Yu)
      cp_X_low_Y_low = X.t() * Y - Xu.t() * Yu; 

      inv_s_delta_Xu = arma::mat( Xu.n_rows, Xu.n_cols);
    } 

    // constructor, precompute Yu, Xu, cp_X_low, cp_X_low_Y_low
    FASTLMM(const arma::vec &Y_, 
            const arma::mat &X_, 
            const arma::vec &Yu_, 
            const arma::mat &Xu_,
            const arma::mat &U_, 
            const arma::vec &s_,
            const arma::mat &cp_X_low_, 
            const arma::mat &cp_X_low_Y_low_){

      Y = Y_;
      X = X_;
      U = U_;
      s = s_;
      Yu = Yu_;
      Xu = Xu_;

      // use weights here

      cp_X_low = cp_X_low_;
      cp_X_low_Y_low = cp_X_low_Y_low_;

      inv_s_delta_Xu = arma::mat( Xu.n_rows, Xu.n_cols);
    } 

    // Accessors
    double get_logLik(){ return logLik; }
    arma::vec get_beta(){ return beta; }
    double get_sigg(){ return sig_g;}
    double get_sige(){ return delta_hat * sig_g;}
    double get_iter(){ return iter;}
    double get_delta(){ return delta_hat;}
    arma::mat get_vcov(){
      return inv(QXX) * sig_g;
    }
    arma::mat get_beta_se(){
      return sqrt(diagvec(get_vcov()));
    }

    // how to combine X and K?
    // to create diagonals of hat matrix?
    double get_edf(); // defined before
    double get_rdf(); // based on Hastie, et al
    arma::vec hatvalues(); // diag of hat matrix
    arma::vec residuals();
    arma::vec predict();
    // df = sum(s[seq_len(rank)]/(s[seq_len(rank)]+delta))

    // compute log likelihood
    double ll(const double &delta);

    void estimate_delta();

    // evaluate logLik, beta, etc at delta value
    void eval_delta( const double &delta){
      logLik = ll( delta );
      delta_hat = delta;
    }

    // Score test
    double score_test(const arma::vec &x_);

    // Update Y, keeping rest constant
    void update_Y( const arma::vec &Y_);

    // Update X, keeping rest constant
    void update_X( const arma::vec &X_);

  private:
    arma::mat Y, Yu;
    arma::mat X, Xu, U;
    arma::vec s;
    arma::mat cp_X_low, cp_X_low_Y_low;
    arma::vec inv_s_delta;
    arma::mat inv_s_delta_Xu;
    arma::mat QXX, QXY;
    arma::mat beta;
    arma::vec r, ru;

    double logLik, sig_g, delta_hat, iter = 0;
};

double FASTLMM::ll(const double &delta ) { 

  double n = X.n_rows;
  double rank = Xu.n_rows;

  inv_s_delta = 1 / (s+delta); 

  // inv_s_delta_Xu   <- inv_s_delta * Xu
  // R recycles over each column,
  //    here do manually
  // inv_s_delta_Xu( Xu.n_rows, Xu.n_cols);
  for(int i=0; i<Xu.n_cols; i++){
    inv_s_delta_Xu.col(i) = inv_s_delta % Xu.col(i);
  }

  // QXX = crossprod(Xu, inv_s_delta_Xu) + cp_X_low / delta
  QXX = Xu.t() * inv_s_delta_Xu + cp_X_low / delta;

  // QXY = crossprod(Xu, inv_s_delta_Yu) + cp_X_low_Y_low / delta
  QXY = Xu.t() * (inv_s_delta % Yu) + cp_X_low_Y_low / delta;

  // beta <<- solve( QXX, QXY)
  beta = arma::solve(QXX, QXY, arma::solve_opts::likely_sympd);

  // # Eval sig_g
  // ru <- Yu - Xu %*% beta
  ru = Yu - Xu * beta;

  // r <- Y - X %*% beta
  r = Y - X * beta;

  // inv_s_delta_ru   <- inv_s_delta * ru
  // vec inv_s_delta_ru = inv_s_delta % ru;

  // Qrr <- crossprod(ru, inv_s_delta_ru) + (crossprod(r)[1] - crossprod(ru)[1])/ delta
  // sig_g <<- Qrr[1] / n
  double QRR = arma::dot(ru, (inv_s_delta % ru)) + (arma::dot(r,r) - arma::dot(ru,ru)) / delta;
  sig_g = QRR / n;

  // use 2.0 to ensure double precision
  double logLik = -n/2.0 * log(2.0*M_PI*sig_g) - 1.0/2.0 * (sum( log(s + delta ) ) + (n-rank) * log(delta)) - n/2.0;

  return logLik;
}
    


// function to be minimized
double ll_alone( double delta_log, void *arg){

  FASTLMM *fit = (FASTLMM *) arg;

  // search is done in log space 
  //  to give faster convergence
  fit->eval_delta( exp(delta_log) );

  return -1.0*fit->get_logLik();
}

void FASTLMM::estimate_delta(){

  double a = -20, b = 20;
  iter = 0;
  
  double max_iter = 100;
  int status;

  // initialize function
  gsl_function F;
  F.function = & ll_alone;
  F.params = this;

  // initialize minimizer
  gsl_min_fminimizer *s;
  s = gsl_min_fminimizer_alloc( gsl_min_fminimizer_brent );

  // Rcout << "a: " << ll_alone(a, this) << std::endl;
  // Rcout << "-4: " << ll_alone(-4, this) << std::endl;
  // Rcout << "b: " << ll_alone(b, this) << std::endl;

  status = gsl_min_fminimizer_set(s, &F, -4, a, b);

  do{
    iter++;
    // Rcout << a << ' ' << b << std::endl;
    status = gsl_min_fminimizer_iterate(s);

    delta_hat = gsl_min_fminimizer_x_minimum(s);
    delta_hat = exp(delta_hat);
    a = gsl_min_fminimizer_x_lower(s);
    b = gsl_min_fminimizer_x_upper(s);

    status = gsl_min_test_interval (a, b, 0.0001, 0.0);
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_min_fminimizer_free(s);
}

#endif
