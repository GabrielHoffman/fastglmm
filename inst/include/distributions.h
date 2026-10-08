/***************************************************************
 * @file    distributions.h
 * @author  Gabriel Hoffman
 * @email   gabriel.hoffman@mssm.edu
 * @brief   Thread-safe distribution functions
 * Copyright (C) 2026 Gabriel Hoffman
 **************************************************************/

#ifndef DISTRIBUTIONS_H_
#define DISTRIBUTIONS_H_

#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/negative_binomial.hpp>

using namespace boost::math;

/* like R::qpois(p, lambda) but threadsafe
*/
// [[Rcpp::export]]
[[maybe_unused]]
static double qpois_boost(const double & p, const double &lambda, const bool &lowertail=true){

  // Initialize Poisson distribution
  boost::math::poisson_distribution<double> dist(lambda);

  // Calculate the quantile (inverse CDF)
  double value = 0;

  if( lowertail ){
    value = quantile(dist, p);
  }else{
    value = quantile(complement(dist, p));
  }

  return value;
}

/* like R::qnbinom(p, size, prob) but threadsafe
*/
// [[Rcpp::export]]
[[maybe_unused]]
static double qnbinom_boost(const double & p, const double &size, const double &prob, const bool &lowertail=true){

  // Initialize NB distribution
  boost::math::negative_binomial_distribution<double> dist(size, prob);

  // Calculate the quantile (inverse CDF)
  double value = 0;

  if( lowertail ){
    value = quantile(dist, p);
  }else{
    value = quantile(complement(dist, p));
  }

  return value;
}


/* like R::dpois(x, lambda) but threadsafe
*/
// [[Rcpp::export]]
[[maybe_unused]]
static double dpois_boost(const double & x, const double &lambda, const bool &useLog=false){

  // Initialize Poisson distribution
  boost::math::poisson_distribution<double> dist(lambda);

  double value = 0;

  if(useLog){
    // Calculate the log PDF
    value = logpdf(dist, x);
  }else{
    // Calculate the PDF
    value = pdf(dist, x);
 }

 return value;
}

/* like R::dnbinom(x, size, prob) but threadsafe
*/
// [[Rcpp::export]]
[[maybe_unused]]
static double dnbinom_boost(const double & x, const double &size, const double &prob, const bool &useLog=false){

  // Initialize NB distribution
  boost::math::negative_binomial_distribution<double> dist(size, prob);

  double value = 0;

  if(useLog){
    // Calculate the log PDF
    value = logpdf(dist, x);
  }else{
    // Calculate the PDF
    value = pdf(dist, x);
 }

 return value;
}


/* Create vector of probabilities for NB distribution

Evaluate NB PDF for y = 0:qmax using recursion to avoid using lgamma()

vec y = generate_seq(qmax);
vec p(y.size());
double prob = theta / (theta + mu);
for(int k = 0; k<y.size(); ++k){
  p(k) = dnbinom_boost(y[k], theta, prob);
}
*/
[[maybe_unused]]
static arma::vec dnbinom_seq_mu_theta(double mu, double theta, int qmax) {
  arma::vec out(qmax + 1);

  double p = mu / (mu + theta);
  double q = theta / (mu + theta);

  // P(Y = 0)
  out[0] = std::pow(q, theta);

  // recurrence:
  // P(y+1) = P(y) * ((y + theta)/(y+1)) * p
  for (int y = 0; y < qmax; ++y) {
    out[y + 1] = out[y] * ( (y + theta) / (y + 1.0) ) * p;
  }

  return out;
}

#endif