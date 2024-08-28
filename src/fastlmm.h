#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// only depends on armadillo 
// Can be imported by plain C++ without Rcpp*
//  by only modifying header above


#ifdef _OPENMP
    // [[Rcpp::plugins(openmp)]]
    #include <omp.h>
#else
    #define omp_get_num_threads() 0
    #define omp_get_thread_num() 0
#endif


#ifndef FASTLMM_H_
#define FASTLMM_H_

using namespace arma;

#include "local_min.h"

namespace fastlmmLib {

class fastlmm_result {       
  public:  
    double logLik, sigSq_g, sigSq_e, delta;
    int iter;
    vec beta, beta_se, weights, ru, r;
    mat vcov;

    fastlmm_result(){}

    fastlmm_result( const double &logLik_,
                    const vec &beta_,
                    const mat &vcov_,
                    const vec &beta_se_,
                    const vec &weights_,
                    const vec &ru_,
                    const vec &r_,
                    const double &delta_,
                    const double &sigSq_g_,
                    const double &sigSq_e_,
                    const int &iter_ ){
      logLik  = logLik_;
      beta    = beta_;
      vcov    = vcov_;
      beta_se = beta_se_;
      weights = weights_;
      ru = ru_;
      r = r_;
      delta   = delta_;
      sigSq_g = sigSq_g_;
      sigSq_e = sigSq_e_;
      iter    = iter_;
    }
};


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
            const vec &weights_);

    // constructor, precompute Yu, Xu
    fastlmm(const T1 &Y_, 
            const T2 &X_, 
            const T3 &U_, 
            const vec &s_,
            const vec &weights_,
            const vec &Yu_, 
            const mat &Xu_);

    // constructor, precompute Yu, Xu, cp_X_low, cp_X_low_Y_low
    fastlmm(const T1 &Y_, 
            const T2 &X_,
            const T3 &U_, 
            const vec &s_,
            const vec &weights_,
            const vec &Yu_, 
            const mat &Xu_,
            const mat &cp_X_low_, 
            const mat &cp_X_low_Y_low_);

    // constructor without response
    fastlmm(const T2 &X_, 
            const T3 &U_, 
            const vec &s_);

    // extract results
    fastlmm_result get_result(){

      mat V = get_vcov();

      vec w = get_weights();

      return fastlmm_result(  get_logLik(),
                              get_beta(),
                              V,
                              sqrt(diagvec(V)),
                              get_weights(),
                              get_ru(),
                              get_r(),
                              get_delta(),
                              get_sigSq_g(),
                              get_sigSq_e(),
                              get_iter());
    }

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
    const double get_rdf(); // based on Hastie, et al
    const vec hatvalues(); // diag of hat matrix
    const vec residuals();
    const vec predict();
    // df = sum(s[seq_len(rank)]/(s[seq_len(rank)]+delta))

    // compute log likelihood
    double ll(const double &delta);

    void estimate_delta(  const double &left,
                          const double &right,
                          const double &tol);

    void update_response(const T1 &Y_, const vec &weights_);
    void update_response(const T1 &Y_,
                         const vec &weights_, 
                         const mat &Yu_);

    vector<fastlmm_result> 
        fit_batch_response(const T1 &Y_all_, 
                           const mat &weights_,
                           const double &delta,
                           const double &left,
                           const double &right,
                           const double &tol,
                           const int &nthreads);

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

  private:
    T1 Y, Yu;
    T2 X, Xu;
    T3 U;
    vec s, weights;
    mat cp_X_low, cp_X_low_Y_low;
    vec inv_s_delta;
    mat inv_s_delta_Xu;
    mat QXX, QXY;
    mat beta;
    vec r, ru;

    double logLik, sigSq_g, delta_hat;
    int iter = 0;
};



// constructor, minimal
template <typename T1, typename T2, typename T3> 
fastlmm<T1, T2, T3>::fastlmm(const T1 &Y_, 
        const T2 &X_, 
        const T3 &U_, 
        const vec &s_,
        const vec &weights_){

  // indicator_decomp
  // modiy this->U  and this->s internally
  // compute sqrt(weights) for 
  // Y <- Y * sqrt(weights)
  // X <- X * sqrt(weights)
  // vec sqrtW = sqrt(weights_);
  // update_weights( Y_, X_, U_, s_, weights_);

  this->Y = Y_;
  this->X = X_;
  this->U = U_;
  this->s = s_;

  // use weights here
  this->weights = weights_;

  this->Yu = U_.t() * Y;
  this->Xu = U_.t() * X;
  this->cp_X_low = X.t() * X - Xu.t() * Xu;
  this->cp_X_low_Y_low = X.t() * Y - Xu.t() * Yu; 
  this->inv_s_delta_Xu = mat( Xu.n_rows, Xu.n_cols);

} 

// constructor, precompute Yu, Xu
template <typename T1, typename T2, typename T3> 
fastlmm<T1, T2, T3>::fastlmm(const T1 &Y_, 
        const T2 &X_, 
        const T3 &U_, 
        const vec &s_,
        const vec &weights_,
        const vec &Yu_, 
        const mat &Xu_){

  this->Y = Y_.t();
  this->X = X_;
  this->U = U_;
  this->s = s_;
  this->weights = weights_;
  this->Yu = Yu_;
  this->Xu = Xu_;
  this->cp_X_low = X.t() * X - Xu.t() * Xu;
  this->cp_X_low_Y_low = X.t() * Y - Xu.t() * Yu; 
  this->inv_s_delta_Xu = mat( Xu.n_rows, Xu.n_cols);

  // use weights here
} 

// constructor, precompute Yu, Xu, cp_X_low, cp_X_low_Y_low
template <typename T1, typename T2, typename T3> 
fastlmm<T1, T2, T3>::fastlmm(const T1 &Y_, 
            const T2 &X_,
            const T3 &U_, 
            const vec &s_,
            const vec &weights_,
            const vec &Yu_, 
            const mat &Xu_,
            const mat &cp_X_low_, 
            const mat &cp_X_low_Y_low_){

  this->Y = Y_.t();
  this->X = X_;
  this->U = U_;
  this->s = s_;
  this->weights = weights_;
  this->Yu = Yu_;
  this->Xu = Xu_;
  this->cp_X_low = cp_X_low_;
  this->cp_X_low_Y_low = cp_X_low_Y_low_;
  this->inv_s_delta_Xu = mat( Xu.n_rows, Xu.n_cols);
  // use weights here?
} 


template <typename T1, typename T2, typename T3> 
fastlmm<T1, T2, T3>::fastlmm( const T2 &X_, 
                  const T3 &U_, 
                  const vec &s_){
  this->X = X_;
  this->U = U_;
  this->s = s_;
  this->Xu = U_.t() * X_;
  this->cp_X_low = X.t() * X - Xu.t() * Xu;
  this->inv_s_delta_Xu = mat( Xu.n_rows, Xu.n_cols);
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
  this->cp_X_low_Y_low = X.t() * Y - Xu.t() * Yu;
}

template <typename T1, typename T2, typename T3>  
double fastlmm<T1, T2, T3>::ll(const double &delta ) { 

  double n = X.n_rows;
  double rank = Xu.n_rows;

  inv_s_delta = 1 / (s+delta);

  // inv_s_delta_Xu   <- inv_s_delta * Xu
  // R recycles over each column,
  //    here do manually
  // inv_s_delta_Xu( Xu.n_rows, Xu.n_cols);
  for(int i=0; i<Xu.n_cols; i++){
    inv_s_delta_Xu.col(i) = inv_s_delta % Xu.col(i);
  };

  // QXX = crossprod(Xu, inv_s_delta_Xu) + cp_X_low / delta
  QXX = Xu.t() * inv_s_delta_Xu + cp_X_low / delta;

  // QXY = crossprod(Xu, inv_s_delta_Yu) + cp_X_low_Y_low / delta
  QXY = Xu.t() * (inv_s_delta % Yu) + cp_X_low_Y_low / delta;

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
double ll_alone_mat( double delta_log, void *arg){

  auto *fit = (fastlmm<mat,mat,mat> *) arg;

  // search is done in log space 
  //  to give faster convergence
  fit->eval_delta( exp(delta_log) );

  return -1.0*fit->get_logLik();
}

// sparse version
double ll_alone_spmat( double delta_log, void *arg){

  auto *fit = (fastlmm<mat,mat,sp_mat> *) arg;

  // search is done in log space 
  //  to give faster convergence
  fit->eval_delta( exp(delta_log) );

  return -1.0*fit->get_logLik();
}

// general case, returns false
template <class T>
bool isSpMatrix(const T &t) { return false;  } 

 // but for sp_mat returns true
template <>
bool isSpMatrix( const sp_mat &t) { return true; } 



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
vector<fastlmm_result> 
  fastlmm<T1, T2, T3>::fit_batch_response( const T1 &Y_all_,
                               const mat &weights_,
                               const double &delta,
                               const double &left,
                               const double &right,
                               const double &tol,
                               const int &nthreads){

    Rcpp::Rcout << "Fit batch response" << std::endl;

  // need to apply weights matrix Y_all_, decomp, and X
  mat Yu_all = U.t() * Y_all_;
  int n_responses = Y_all_.n_cols;

  // store results
  vector<fastlmm_result> result(n_responses, fastlmm_result());

  // NOTE: Do not use Rcpp in parallel section
  // "C stack usage is too close to the limit"

  #ifdef _OPENMP 
  // set threads
  omp_set_num_threads(nthreads);
  // disable nested parallelism
  omp_set_max_active_levels(1);
  #endif

  #pragma omp parallel
  {
    // initialize
    fastlmm fit = fastlmm(X, U, s);

    // iterate through responses 
    #pragma omp for 
    for( int i = 0; i < n_responses; i++){

      fit.update_response(Y_all_.col(i), 
                          weights_.col(i), 
                          Yu_all.col(i));

      if( delta > 0 ){
        fit.eval_delta( delta ); 
      }else{
        fit.estimate_delta( left, right, tol );
      }

      #pragma omp critical
      result.at(i) = fit.get_result();
    }
  }

  return result;
}

}


#endif
