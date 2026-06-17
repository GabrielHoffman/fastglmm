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
double qpois_boost(const double & p, const double &lambda, const bool &lowertail=true){

  // Initialize Poisson distribution
  boost::math::poisson_distribution<double> dist(lambda);

  // Calculate the quantile (inverse CDF)
  double value;

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
double qnbinom_boost(const double & p, const double &size, const double &prob, const bool &lowertail=true){

  // Initialize NB distribution
  boost::math::negative_binomial_distribution<double> dist(size, prob);

  // Calculate the quantile (inverse CDF)
  double value;

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
double dpois_boost(const double & x, const double &lambda){

  // Initialize Poisson distribution
  boost::math::poisson_distribution<double> dist(lambda);

  // Calculate the PDF
  return pdf(dist, x);
}

/* like R::dnbinom(x, size, prob) but threadsafe
*/
// [[Rcpp::export]]
double dnbinom_boost(const double & x, const double &size, const double &prob){

  // Initialize NB distribution
  boost::math::negative_binomial_distribution<double> dist(size, prob);

  // Calculate the PDF
  return pdf(dist, x);
}


#endif