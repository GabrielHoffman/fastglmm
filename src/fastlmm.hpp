#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// only depends on armadillo and gsl
// Can be imported by plain C++ without Rcpp*
//  by only modifying header above

#ifndef FASTLMM_H_
#define FASTLMM_H_

#include <RcppGSL.h>

#include <gsl/gsl_min.h>
#include <gsl/gsl_errno.h>


class FASTLMM_result {       
  public:  
    double logLik, sig_g, sig_e, delta;
    int iter;
    arma::vec beta, beta_se;
    arma::mat vcov;

    FASTLMM_result(){}

    FASTLMM_result( const double &logLik_,
                    const arma::vec &beta_,
                    const arma::mat &vcov_,
                    const arma::vec &beta_se_,
                    const double &delta_,
                    const double &sig_g_,
                    const double &sig_e_,
                    const int &iter_ ){
      logLik = logLik_;
      beta   = beta_;
      vcov   = vcov_;
      beta_se= beta_se_;
      delta  = delta_;
      sig_g  = sig_g_;
      sig_e  = sig_e_;
      iter   = iter_;
    }
};

template <typename T> 
class FASTLMM {       
  public:  

    // constructor, minimal
    // template <typename T> 
    FASTLMM(const arma::vec &Y_, 
            const arma::mat &X_, 
            const T &U_, 
            const arma::vec &s_);

    // constructor, precompute Yu, Xu
    // template <typename T> 
    FASTLMM(const arma::vec &Y_, 
            const arma::mat &X_, 
            const T &U_, 
            const arma::vec &s_,
            const arma::vec &Yu_, 
            const arma::mat &Xu_);

    // constructor, precompute Yu, Xu, cp_X_low, cp_X_low_Y_low
    // template <typename T> 
    FASTLMM(const arma::vec &Y_, 
            const arma::mat &X_,
            const T &U_, 
            const arma::vec &s_,
            const arma::vec &Yu_, 
            const arma::mat &Xu_,
            const arma::mat &cp_X_low_, 
            const arma::mat &cp_X_low_Y_low_);

    // constructor with out response
    // template <typename T> 
    FASTLMM(const arma::mat &X_, 
            const T &U_, 
            const arma::vec &s_);

    // extract results
    FASTLMM_result get_result(){

      arma::mat V = get_vcov();

      return FASTLMM_result(  get_logLik(),
                              get_beta(),
                              V,
                              sqrt(diagvec(V)),
                              get_delta(),
                              get_sigg(),
                              get_sige(),
                              get_iter());
    }

    // Accessors
    const double get_logLik(){ return this->logLik; }
    const arma::vec get_beta(){ return this->beta; }
    const double get_sigg(){ return sig_g;}
    const double get_sige(){ 
      return this->delta_hat * this->sig_g;
    }
    const double get_iter(){ return this->iter;}
    const double get_delta(){ return this->delta_hat;}
    const arma::mat get_vcov(){
      return inv(this->QXX) * this->sig_g;
    }
    const arma::mat get_beta_se(){
      return sqrt(diagvec(get_vcov()));
    }

    // how to combine X and K?
    // to create diagonals of hat matrix?
    const double get_edf(); // defined before
    const double get_rdf(); // based on Hastie, et al
    const arma::vec hatvalues(); // diag of hat matrix
    const arma::vec residuals();
    const arma::vec predict();
    // df = sum(s[seq_len(rank)]/(s[seq_len(rank)]+delta))

    // compute log likelihood
    double ll(const double &delta);

    void estimate_delta();

    void update_response(const arma::vec &Y_);
    void update_response(const arma::vec &Y_, 
                         const arma::vec &Yu_);

    std::vector<FASTLMM_result> 
        fit_batch_response(const arma::mat &Y_all_, 
                           const double &delta);

    // evaluate logLik, beta, etc at delta value
    void eval_delta( const double &delta){
      this->logLik = ll( delta );
      this->delta_hat = delta;
    }

    // Score test
    double score_test(const arma::vec &x_);

    // Update Y, keeping rest constant
    void update_Y( const arma::vec &Y_);

    // Update X, keeping rest constant
    void update_X( const arma::vec &X_);

  private:
    arma::vec Y, Yu;
    arma::mat X, Xu;
    T U;
    arma::vec s;
    arma::mat cp_X_low, cp_X_low_Y_low;
    arma::vec inv_s_delta;
    arma::mat inv_s_delta_Xu;
    arma::mat QXX, QXY;
    arma::mat beta;
    arma::vec r, ru;

    double logLik, sig_g, delta_hat, iter = 0;
};



// constructor, minimal
template <typename T> 
FASTLMM<T>::FASTLMM(const arma::vec &Y_, 
        const arma::mat &X_, 
        const T &U_, 
        const arma::vec &s_){

  this->Y = Y_;
  this->X = X_;
  this->U = U_;
  this->s = s_;
  this->Yu = U_.t() * Y_;
  this->Xu =  U_.t() * X_;
  this->cp_X_low = X.t() * X - Xu.t() * Xu;
  this->cp_X_low_Y_low = X.t() * Y - Xu.t() * Yu; 
  this->inv_s_delta_Xu = arma::mat( Xu.n_rows, Xu.n_cols);
  // use weights here
} 

// constructor, precompute Yu, Xu
template <typename T> 
FASTLMM<T>::FASTLMM(const arma::vec &Y_, 
        const arma::mat &X_, 
        const T &U_, 
        const arma::vec &s_,
        const arma::vec &Yu_, 
        const arma::mat &Xu_){

  this->Y = Y_;
  this->X = X_;
  this->U = U_;
  this->s = s_;
  this->Yu = Yu_;
  this->Xu = Xu_;
  this->cp_X_low = X.t() * X - Xu.t() * Xu;
  this->cp_X_low_Y_low = X.t() * Y - Xu.t() * Yu; 
  this->inv_s_delta_Xu = arma::mat( Xu.n_rows, Xu.n_cols);

  // use weights here
} 

// constructor, precompute Yu, Xu, cp_X_low, cp_X_low_Y_low
template <typename T> 
FASTLMM<T>::FASTLMM(const arma::vec &Y_, 
            const arma::mat &X_,
            const T &U_, 
            const arma::vec &s_,
            const arma::vec &Yu_, 
            const arma::mat &Xu_,
            const arma::mat &cp_X_low_, 
            const arma::mat &cp_X_low_Y_low_){

  this->Y = Y_;
  this->X = X_;
  this->U = U_;
  this->s = s_;
  this->Yu = Yu_;
  this->Xu = Xu_;
  this->cp_X_low = cp_X_low_;
  this->cp_X_low_Y_low = cp_X_low_Y_low_;
  this->inv_s_delta_Xu = arma::mat( Xu.n_rows, Xu.n_cols);
  // use weights here?
} 

template <typename T> 
FASTLMM<T>::FASTLMM( const arma::mat &X_, 
                  const T &U_, 
                  const arma::vec &s_){
  this->X = X_;
  this->U = U_;
  this->s = s_;
  this->Xu = U_.t() * X_;
  this->cp_X_low = X.t() * X - Xu.t() * Xu;
  this->inv_s_delta_Xu = arma::mat( Xu.n_rows, Xu.n_cols);
}

template <typename T> 
void FASTLMM<T>::update_response(const arma::vec &Y_){

  update_response(Y, U.t() * Y_);
} 

template <typename T> 
void FASTLMM<T>::update_response(const arma::vec &Y_, 
                              const arma::vec &Yu_){
  this->Y = Y_;
  this->Yu = Yu_;  
  this->cp_X_low_Y_low = X.t() * Y - Xu.t() * Yu;
}


template <typename T> 
double FASTLMM<T>::ll(const double &delta ) { 

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

  // Qrr <- crossprod(ru, inv_s_delta_ru) + (crossprod(r)[1] - crossprod(ru)[1])/ delta
  // sig_g <<- Qrr[1] / n
  double QRR = arma::dot(ru, (inv_s_delta % ru)) + (arma::dot(r,r) - arma::dot(ru,ru)) / delta;
  sig_g = QRR / n;

  // use 2.0 to ensure double precision
  double logLik = -n/2.0 * log(2.0*M_PI*sig_g) - 1.0/2.0 * (sum( log(s + delta ) ) + (n-rank) * log(delta)) - n/2.0;

  return logLik;
}
    


// function to be minimized
double ll_alone_mat( double delta_log, void *arg){

  FASTLMM<arma::mat> *fit = (FASTLMM<arma::mat> *) arg;

  // search is done in log space 
  //  to give faster convergence
  fit->eval_delta( exp(delta_log) );

  return -1.0*fit->get_logLik();
}
// sparse version
double ll_alone_spmat( double delta_log, void *arg){

  FASTLMM<arma::sp_mat> *fit = (FASTLMM<arma::sp_mat> *) arg;

  // search is done in log space 
  //  to give faster convergence
  fit->eval_delta( exp(delta_log) );

  return -1.0*fit->get_logLik();
}

// general case, returns false
template <class T>
bool isSpMatrix(const T &t) { return false;  } 

 // but for arma::sp_mat returns true
template <>
bool isSpMatrix( const arma::sp_mat &t) { return true; } 


template <typename T>
void FASTLMM<T>::estimate_delta(){

  double a = -20, b = 20;
  iter = 0;
  
  double max_iter = 100;
  int status;

  // initialize function
  gsl_function F;  
  F.params = this;

  // Since F.function can't take templated function
  if( isSpMatrix( U ) ){
    F.function = & ll_alone_spmat;
  }else{
    F.function = & ll_alone_mat;
  }

  // initialize minimizer
  gsl_min_fminimizer *s;
  s = gsl_min_fminimizer_alloc( gsl_min_fminimizer_brent );

  // Rcout << "a: " << ll_alone(a, this) << std::endl;
  // Rcout << "-4: " << ll_alone(-4, this) << std::endl;
  // Rcout << "b: " << ll_alone(b, this) << std::endl;

  // if this fails, itertively half initial value
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




template <typename T>
std::vector<FASTLMM_result> 
  FASTLMM<T>::fit_batch_response( const arma::mat &Y_all_,
                               const double &delta){

  // responses are stored as _rows_
  arma::mat Yu_all = Y_all_ * U;

  int n_responses = Y_all_.n_rows;

  // store results
  std::vector<FASTLMM_result> result(n_responses, FASTLMM_result()); 
  int OMP_CHUNK_SIZE = n_responses / omp_get_num_threads();

  // Rcout << "omp_get_num_threads: " << omp_get_num_threads() << std::endl;
  // Rcout << "OMP_CHUNK_SIZE: " << OMP_CHUNK_SIZE << std::endl;

  // disable nested parallelism
  omp_set_nested(0);
  #pragma omp parallel
  {
    // initialize
    FASTLMM fit = FASTLMM(X, U, s);

    // iterate thru responses i.e. rows
    #pragma omp for schedule(static, OMP_CHUNK_SIZE)
    for( int i = 0; i < n_responses; i++){
      fit.update_response(Y_all_.row(i).t(), Yu_all.row(i).t());

      if( delta > 0 ){
        fit.eval_delta( delta ); 
      }else{
        fit.estimate_delta();
      }

      #pragma omp critical
      result.at(i) = fit.get_result();
    }
  }

  return result;
}


// List fastlmm_batchX_c()

#endif
