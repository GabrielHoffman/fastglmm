#ifndef LOG_MOMENTS_NB_H_
#define LOG_MOMENTS_NB_H_

// if -D USE_R, use RcppArmadillo library
#ifdef USE_R
// [[Rcpp::depends(RcppParallel)]]  
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#else
#include <armadillo>
#include <tbb/tbb.h>
#endif

#include <Rcpp.h>
#include <utility> 

#include <distributions.h> 

using namespace arma;
using namespace std;

vec generate_seq(const int &qmax) {
  // Initialize an IntegerVector of size (qmax + 1)
  vec y(qmax + 1);
  
  // Fill the vector with values starting from 0 to qmax
  std::iota(y.begin(), y.end(), 0);
  
  return y;
}

std::pair<double, double> _log_moments_nb_mu(
                    const vec & mu,
                    const double &theta, // overdispersion
                    const double &c = 1.0, // pseudocount
                    const double &p_tail = 1e-4) {

  int n = mu.size();

  vec mean_out(n);
  vec var_out(n);

  // for each entry in mu
  for(int i = 0; i < n; ++i) {

    double qmax_d;

    // get quantiles
    // can be off by 1 due to differences in rounding
    // but the boost functions are thread-safe
    if( isinf(theta) ){
      // qmax_d = R::qpois(p_tail, mu[i], false, false);
      qmax_d = qpois_boost(p_tail, mu[i], false);
    }else{
      // qmax_d = R::qnbinom(p_tail, theta, theta / (theta + mu[i]), false, false);      
      qmax_d = qnbinom_boost(p_tail, theta, theta / (theta + mu[i]), false);
    }

    int qmax = static_cast<int>(qmax_d);

    // y <- seq(0, qmax)
    vec y = generate_seq(qmax);

    vec p(y.size());

    if( isinf(theta) ){      

      // Poisson
      for(int k = 0; k<y.size(); k++){
        p(k) = dpois_boost(y[k], mu[i]);
      }
    }else{
      // Negative binomial
      double prob = theta / (theta + mu[i]);
      for(int k = 0; k<y.size(); k++){
       p(k) = dnbinom_boost(y[k], theta, prob);
      }
    }

    // Normalize probabilies to sum to 1
    p = p / sum(p);

    // Compute log response plus pseudocount
    vec z = log(y + c);

    // Expected value of z
    double ez = sum(z % p);

    // Variance of z
    // vz <- sum((z - ez)^2 * p)
    vec dz = z - ez;
    double vz = sum(dz % dz % p);

    mean_out[i] = ez;
    var_out[i] = vz;
  }

  double signal = var(mean_out);
  double noise = mean(var_out);

  return {signal, noise};
}



/* Log Moments of NB given X, Beta, and offset */
tuple<vec, vec, vec> _log_moments_nb_XB(
  const mat & X,
  const mat & Beta,
  const vec & offset,
  const vec & theta, // overdispersion
  const string &method,
  const double & c = 1.0, // pseudocount
  const double & p_tail = 1e-4,
  const int & nthreads = 10) {

  int n_responses = Beta.n_rows;
  vec signal(n_responses);
  vec noise(n_responses);
  vec alpha(n_responses);

  // Parallel part using Thread Building Blocks
  tbb::task_arena limited_arena(nthreads);
  limited_arena.execute([&] {
  tbb::parallel_for(
    tbb::blocked_range<int>(0, n_responses, 10), 
    [&](const tbb::blocked_range<int>& r){ 

    disable_parallel_blas();

    vec mu(X.n_rows);
    vec var_poisson(X.n_rows), var_overdisp(X.n_rows);
    for (int j = r.begin(); j != r.end(); ++j) {

      // create mu
      mu = exp(X * Beta.row(j).t() + offset);

      // Poisson variance
      var_poisson = mu / square(mu+c);

      // Overdispersion variance
      var_overdisp = (square(mu)/theta(j)) / square(mu+c);

      // Poisson noise ratio
      alpha(j) = mean( var_poisson / (var_poisson + var_overdisp));

      if( method == "exact" && mean(mu) < 1e8){
        // evaluate moments
        auto [signal_, noise_] = _log_moments_nb_mu(mu, theta(j), c, p_tail);    
        signal(j) = signal_;
        noise(j) = noise_; 
      }else{        
        signal(j) = var(log(mu + c));
        noise(j) = mean(var_poisson + var_overdisp); 
      }
    }
   }); });

  return {signal, noise, alpha};
}


/* Log Moments of NB given BLUP, X, Beta, Z, offset */
template <typename T> 
tuple<vec, vec, vec> _log_moments_nb_BlupXBZ(
  const mat & BLUP,
  const mat & X,
  const mat & Beta,
  const T & Z,  
  const vec & weights,
  const vec & offset,
  const vec & theta, // overdispersion
  const vec & delta, // variance component ratio
  const string &method,
  const string &dcmpMethod, 
  const double & c = 1.0, // pseudocount
  const double & p_tail = 1e-4,
  const int & nthreads = 10) {

  int n_responses = Beta.n_rows;
  vec signal(n_responses);
  vec noise(n_responses);
  vec alpha(n_responses);

  uvec idx_drop = find(weights == 0.0);
  double n_active = weights.n_elem - idx_drop.n_elem;

  ZTYPE type = dcmpMethod == "categorical" ? 
                CATEGORICAL : GENERAL;
  spectralDecomp<T> dcmp(Z, type);

  // Parallel part using Thread Building Blocks
  tbb::task_arena limited_arena(nthreads);
  limited_arena.execute([&] {
  tbb::parallel_for(
    tbb::blocked_range<int>(0, n_responses, 10), 
    [&](const tbb::blocked_range<int>& r){ 

    disable_parallel_blas();

    // need local version due to reweighting
    spectralDecomp<T> dcmp_lcl(dcmp);

    vec mu(X.n_rows), eta(X.n_rows);
    vec var_poisson(X.n_rows), var_overdisp(X.n_rows);
    T A;

    for (int j = r.begin(); j != r.end(); ++j) { 

      shared_ptr<GLMFamily> fam = getGLMFamily( "nb:" + to_string(theta[j]) );

      // create eta
      eta = offset + X * Beta.row(j).t() + Z * BLUP.col(j);
      mu = exp(eta);

      // Poisson variance
      var_poisson = mu / square(mu+c);

      // Overdispersion variance
      var_overdisp = (square(mu)/theta(j)) / square(mu+c);

      // Poisson ratio
      alpha(j) = mean( var_poisson / (var_poisson + var_overdisp));

      if( method == "exact" && mean(mu) < 1e8){
        // evaluate moments
        auto [signal_, noise_] = _log_moments_nb_mu(mu, theta(j), c, p_tail);    
        signal(j) = signal_;
        noise(j) = noise_; 
      }else{        
        signal(j) = var(log(mu + c));
        noise(j) = mean(var_poisson + var_overdisp); 
      }
    }
   }); });

  return {signal, noise, alpha};
}





#endif








