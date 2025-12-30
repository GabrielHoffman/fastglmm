#pragma once
#include <cmath>
#include <limits>
#include <algorithm>

// AS241: inverse standard normal CDF (quantile).
// Returns z such that P(Z <= z) = p  (if lower_tail=true).
// Supports log_p (p given as log(p)).
inline double qnorm_as241(double p,
                          double mean = 0.0,
                          double sd = 1.0,
                          bool lower_tail = true,
                          bool log_p = false)
{
  if (!(sd > 0.0) || !std::isfinite(sd)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  // Handle log_p and tail choice
  double pp;
  if (log_p) {
    if (p > 0.0) return std::numeric_limits<double>::quiet_NaN(); // log(p) must be <= 0
    pp = std::exp(p);
  } else {
    pp = p;
  }

  if (!lower_tail) pp = 1.0 - pp;

  // Bounds
  if (pp <= 0.0) return -std::numeric_limits<double>::infinity();
  if (pp >= 1.0) return  std::numeric_limits<double>::infinity();

  // Coefficients (AS241)
  static const double a[8] = {
    3.3871328727963666080,
    1.3314166789178437745e+2,
    1.9715909503065514427e+3,
    1.3731693765509461125e+4,
    4.5921953931549871457e+4,
    6.7265770927008700853e+4,
    3.3430575583588128105e+4,
    2.5090809287301226727e+3
  };
  static const double b[8] = {
    1.0,
    4.2313330701600911252e+1,
    6.8718700749205790830e+2,
    5.3941960214247511077e+3,
    2.1213794301586595867e+4,
    3.9307895800092710610e+4,
    2.8729085735721942674e+4,
    5.2264952788528545610e+3
  };
  static const double c[8] = {
    1.42343711074968357734,
    4.63033784615654529590,
    5.76949722146069140550,
    3.64784832476320460504,
    1.27045825245236838258,
    2.41780725177450611770e-1,
    2.27238449892691845833e-2,
    7.74545014278341407640e-4
  };
  static const double d[8] = {
    1.0,
    2.05319162663775882187,
    1.67638483018380384940,
    6.89767334985100004550e-1,
    1.48103976427480074590e-1,
    1.51986665636164571966e-2,
    5.47593808499534494600e-4,
    1.05075007164441684324e-9
  };
  static const double e[8] = {
    6.65790464350110377720,
    5.46378491116411436990,
    1.78482653991729133580,
    2.96560571828504891230e-1,
    2.65321895265761230930e-2,
    1.24266094738807843860e-3,
    2.71155556874348757815e-5,
    2.01033439929228813265e-7
  };
  static const double f[8] = {
    1.0,
    5.99832206555887937690e-1,
    1.36929880922735805310e-1,
    1.48753612908506148525e-2,
    7.86869131145613259100e-4,
    1.84631831751005468180e-5,
    1.42151175831644588870e-7,
    2.04426310338993978564e-15
  };

  auto poly = [](const double* coef, int n, double x) {
    // Horner: coef[0] + coef[1] x + ... + coef[n-1] x^(n-1)
    double v = coef[n - 1];
    for (int i = n - 2; i >= 0; --i) v = v * x + coef[i];
    return v;
  };

  constexpr double split1 = 0.425;
  constexpr double split2 = 5.0;
  constexpr double const1 = 0.180625;
  constexpr double const2 = 1.6;

  double q = pp - 0.5;
  double z;

  if (std::fabs(q) <= split1) {
    // Central region
    double r = const1 - q*q;
    z = q * poly(a, 8, r) / poly(b, 8, r);
  } else {
    // Tails
    double r = (q < 0.0) ? pp : (1.0 - pp);
    r = std::sqrt(-std::log(r));

    if (r <= split2) {
      r -= const2;
      z = poly(c, 8, r) / poly(d, 8, r);
    } else {
      r -= split2;
      z = poly(e, 8, r) / poly(f, 8, r);
    }
    if (q < 0.0) z = -z;
  }

  return mean + sd * z;
}
