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

// for values greater than this, use approximation
#define COUNT_CEILING 1e8

// Evaluate poisson pdf across many values of y using Panjer recursion
std::pair<double, double> _log_moments_poisson_exact_fast(
      const double& m,
      const double &c = 1.0,
      const double &p_tail = 1e-4){

  // For large count values, use approximation
  if( m > COUNT_CEILING ){
    double vz = boost::math::trigamma(m);
    double ez = log(m); 

    return {ez, vz};
  }

  int qmax = qpois_boost(p_tail, m, false);
  if (qmax < 1) qmax = 1;

  int mode = std::floor(m);
  if (mode > qmax) mode = qmax;

  double sum_p  = 0.0;
  double sum_z  = 0.0;
  double sum_z2 = 0.0;

  // anchor at mode with arbitrary scale
  double py = 1.0;

  // include mode
  if (!(c == 0.0 && mode == 0)) {
    double z = std::log(mode + c);
    sum_p  += py;
    sum_z  += py * z;
    sum_z2 += py * z * z;
  }

  // recurse downward from mode to 0:
  // p_{y-1} = p_y * y / mu
  py = 1.0;
  for (int y = mode; y >= 1; --y) {
    py *= y / m;  // now py is scaled p_{y-1}

    int yy = y - 1;

    if (c == 0.0 && yy == 0) continue;

    double z = std::log(yy + c);
    sum_p  += py;
    sum_z  += py * z;
    sum_z2 += py * z * z;
  }

  // recurse upward from mode to qmax:
  // p_y = p_{y-1} * mu / y
  py = 1.0;
  for (int y = mode + 1; y <= qmax; ++y) {
    py *= m / y;  // now py is scaled p_y

    double z = std::log(y + c);
    sum_p  += py;
    sum_z  += py * z;
    sum_z2 += py * z * z;
  }

  double ez = sum_z / sum_p;
  double vz = sum_z2 / sum_p - ez * ez;

  return {ez, vz};
}

std::pair<vec, vec> _log_moments_poisson_exact_fast(
      const arma::vec& mu,
      const double &c = 1.0,
      const double &p_tail = 1e-4) {

  int n = mu.n_elem;
  arma::vec mean_out(n);
  arma::vec var_out(n);

  for (int i = 0; i < n; ++i) {

    auto [ez, vz] = _log_moments_poisson_exact_fast(mu[i], c, p_tail);
    
    mean_out[i] = ez;
    var_out[i] = vz;
  }

  // signal = var(mean_out);
  // noise = mean(var_out);
  return {mean_out, var_out};
}


// Evaluate NB PDF for y = 0:qmax using recursion Panjer recursion
std::pair<double, double> _log_moments_nb_exact_fast(
      const double& m,
      const double &theta,
      const double &c = 1.0,
      const double &p_tail = 1e-4) {

  // For large count values, use approximation
  if( m > COUNT_CEILING ){
    double vz = boost::math::trigamma(1/(1/m + 1/theta));
    double ez = log(m + c) - 0.5*(m + pow(m,2)/theta)/pow(m+c,2);

    return {ez, vz};
  }

  // R parameterization: size = theta, mu = m
  double prob = theta / (theta + m);

  // Use R's qnbinom for truncation cutoff
  // int qmax = R::qnbinom(p_tail, theta, prob, false, false);
  int qmax = qnbinom_boost(p_tail, theta, prob, false);
  if (qmax < 1) qmax = 1;

  double p = m / (m + theta);

   // NB mode: floor((theta - 1) * p / (1 - p)), for theta > 1
  // In mean/size form this is floor((theta - 1) * mu / theta).
  int mode = 0;
  if (theta > 1.0) {
    mode = std::floor((theta - 1.0) * m / theta);
    if (mode < 0) mode = 0;
    if (mode > qmax) mode = qmax;
  }

  double sum_p  = 0.0;
  double sum_z  = 0.0;
  double sum_z2 = 0.0;

  // Anchor at mode with arbitrary scaled mass
  double py = 1.0;

  // Include mode
  if (!(c == 0.0 && mode == 0)) {
    double z = std::log(mode + c);
    sum_p  += py;
    sum_z  += py * z;
    sum_z2 += py * z * z;
  }

  // Downward recurrence:
  // p_{y-1} = p_y * y / ((y + theta - 1) * p)
  py = 1.0;
  for (int y = mode; y >= 1; --y) {
    py *= y / ((y + theta - 1.0) * p);

    int yy = y - 1;

    if (c == 0.0 && yy == 0) continue;

    double z = std::log(yy + c);
    sum_p  += py;
    sum_z  += py * z;
    sum_z2 += py * z * z;
  }

  // Upward recurrence:
  // p_y = p_{y-1} * ((y + theta - 1) / y) * p
  py = 1.0;
  for (int y = mode + 1; y <= qmax; ++y) {
    py *= ((y + theta - 1.0) / y) * p;

    double z = std::log(y + c);
    sum_p  += py;
    sum_z  += py * z;
    sum_z2 += py * z * z;
  }

  double ez = sum_z / sum_p;
  double vz = sum_z2 / sum_p - ez * ez;

  return {ez, vz};
}


std::tuple<vec, vec, double> _log_moments_nb_exact_fast(
      const arma::vec& mu,
      const double &theta,
      const double &c = 1.0,
      const double &p_tail = 1e-4) {

  // Evaluate Poisson signal and noise
  auto [ez_poisson, vz_poisson] = _log_moments_poisson_exact_fast(mu, c, p_tail);
  
  // Evaluate NB signal and noise
  int n = mu.n_elem;
  arma::vec mean_out(n);
  arma::vec var_out(n);

  for (int i = 0; i < n; ++i) {
   
    auto [ez, vz] = _log_moments_nb_exact_fast(mu[i], theta, c, p_tail);

    mean_out[i] = ez;
    var_out[i] = vz;
  }

  // Poisson noise ratio
  double alpha = mean( vz_poisson / var_out );

  // cap alpha at 1
  // empirical fraction can exceed 1 due to delta approx
  alpha = std::min(alpha, 1.0);

  return {mean_out, var_out, alpha};
}

std::tuple<double, double, double> _log_moments_nb_mu(
                    const vec & mu,
                    const double &theta, // overdispersion
                    const string &method,
                    const double &c = 1.0, // pseudocount
                    const double &p_tail = 1e-4) {

  double signal, noise, alpha;

  if( method == "approximate"){
    vec sq_mu_c = square(mu+c);
    vec sq_mu_theta = square(mu)/theta;

    // Poisson variance
    // vec var_poisson = mu / square(mu+c);
    vec var_poisson = mu / sq_mu_c;

    // Overdispersion variance
    // vec var_overdisp = (square(mu)/theta) / square(mu+c);
    vec var_overdisp = sq_mu_theta / sq_mu_c;

    // Poisson noise ratio
    vec vz = var_poisson + var_overdisp;
    alpha = mean( var_poisson / vz);

    // vec ez = log(mu + c) - 0.5*(mu+square(mu)/theta)/square(mu+c);
    vec ez = log(mu + c) - 0.5*(mu + sq_mu_theta)/sq_mu_c;

    signal = var(ez);
    noise = mean(vz);
  }else{
    // Exact
    if( isinf(theta) ){
      // Poisson
      auto [ez, vz] = _log_moments_poisson_exact_fast(mu, c, p_tail);

      signal = var(ez);
      noise = mean(vz);
      alpha = 1.0;
    }else{
      // Negative binomial
      auto [ez, vz, a] = _log_moments_nb_exact_fast(mu, theta, c, p_tail);

      signal = var(ez);
      noise = mean(vz);
      alpha = a;
    }
  }

  return {signal, noise, alpha};
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
    for (int j = r.begin(); j != r.end(); ++j) {

      // create mu
      mu = exp(X * Beta.row(j).t() + offset);

      // evaluate moments
      auto [_signal, _noise, _alpha] = _log_moments_nb_mu(mu, theta(j), method, c, p_tail);    
  
      signal(j) = _signal;
      noise(j) = _noise; 
      alpha(j) = _alpha; 
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

    for (int j = r.begin(); j != r.end(); ++j) { 

      shared_ptr<GLMFamily> fam = getGLMFamily( "nb:" + to_string(theta[j]) );

      // create eta
      eta = offset + X * Beta.row(j).t() + Z * BLUP.col(j);
      mu = exp(eta);

      // evaluate moments
      auto [_signal, _noise, _alpha] = _log_moments_nb_mu(mu, theta(j), method, c, p_tail);    
  
      signal(j) = _signal;
      noise(j) = _noise; 
      alpha(j) = _alpha; 
    }
   }); });

  return {signal, noise, alpha};
}




/*
vec generate_seq(const int &qmax) {
  // Initialize an IntegerVector of size (qmax + 1)
  vec y(qmax + 1);
  
  // Fill the vector with values starting from 0 to qmax
  std::iota(y.begin(), y.end(), 0);
  
  return y;
}


std::pair<double, double> _log_moments_nb_exact(
  const double &mu, 
  const double &theta,
  const double &c = 1.0, // pseudocount
  const double &p_tail = 1e-4){

  double qmax_d;

  // get quantiles
  // can be off by 1 due to differences in rounding
  // but the boost functions are thread-safe
  if( isinf(theta) ){
    // qmax_d = R::qpois(p_tail, mu, false, false);
    qmax_d = qpois_boost(p_tail, mu, false);
  }else{
    // qmax_d = R::qnbinom(p_tail, theta, theta / (theta + mu), false, false);      
    qmax_d = qnbinom_boost(p_tail, theta, theta / (theta + mu), false);
  }

  int qmax = static_cast<int>(qmax_d);

  vec p = dnbinom_seq_mu_theta(mu, theta, qmax);

  // y <- seq(0, qmax)
  vec y = generate_seq(qmax);

  if( c == 0){    
    // if pseudocount is zero, drop y==0 element
    y.shed_row(0);
    p.shed_row(0);
  }

  // Normalize probabilies to sum to 1
  p = p / sum(p);

  // Compute log response plus pseudocount
  vec z = log(y + c);

  // Expected value of z
  double ez = sum(z % p);

  // Variance of z
  vec dz = z - ez;
  double vz = sum(dz % dz % p);

  return {ez, vz};
}*/

/*
std::pair<double, double> _log_moments_nb_exact(
  const double &mu, 
  const double &theta,
  const double &c = 1.0, // pseudocount
  const double &p_tail = 1e-4){

  double qmax_d;

  // get quantiles
  // can be off by 1 due to differences in rounding
  // but the boost functions are thread-safe
  if( isinf(theta) ){
    // qmax_d = R::qpois(p_tail, mu, false, false);
    qmax_d = qpois_boost(p_tail, mu, false);
  }else{
    // qmax_d = R::qnbinom(p_tail, theta, theta / (theta + mu), false, false);      
    qmax_d = qnbinom_boost(p_tail, theta, theta / (theta + mu), false);
  }

  int qmax = static_cast<int>(qmax_d);

  // y <- seq(0, qmax)
  vec y = generate_seq(qmax);

  if( c == 0){    
    // if pseudocount is zero, drop y==0 element
    y.shed_row(0);
  }

  vec p(y.size());
  if( isinf(theta) ){  
    // Poisson
    for(int k = 0; k<y.size(); k++){
      p(k) = dpois_boost(y[k], mu);
    }
  }else{
    // Negative binomial
    double prob = theta / (theta + mu);

    for(int k = 0; k<y.size(); ++k){
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
  vec dz = z - ez;
  double vz = sum(dz % dz % p);

  return {ez, vz};
}*/


// // [[Rcpp::export]]
// int f( const vec & mu,
//                     const double &theta, // overdispersion
//                     const double &c = 1.0, // pseudocount
//                     const double &p_tail = 1e-4) {

//   for(int i = 0; i < mu.size(); ++i) {

//     pair<double, double> res = _log_moments_nb_exact(mu[i], theta, c, p_tail);    
//   }
//   return 1;
// }


#endif








