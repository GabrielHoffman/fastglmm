# Estimate theta of the Negative Binomial

Given the estimated mean vector, estimate `theta` overdispersion
parameter of the negative binomial distribution.

## Usage

``` r
nb_theta(
  y,
  mu,
  n,
  weights,
  left = -10,
  right = 20,
  tol = .Machine$double.eps^0.25
)
```

## Arguments

- y:

  vector of observed values from the negative binomial.

- mu:

  estimated mean vector from generalized linear model

- n:

  number of data points (defaults to the sum of `weights`)

- weights:

  sample weights. If missing, set to 1

- left:

  left boundary for estimating theta on a log scale. So -10 corresponds
  to a minimum theta value of exp(-10)

- right:

  right boundary for estimating theta on a log scale

- tol:

  tolerance of optimization

## Details

Estimate overdispersion parameter for negative binomial distribution
using univariate optimization with Brent's method. For very large
datasets, this can be much faster than the Newton-Raphson method used by
[`MASS::theta.ml()`](https://rdrr.io/pkg/MASS/man/theta.md.html). This
is faster since the [`lgamma()`](https://rdrr.io/r/base/Special.html)
function used in the likelihood is faster than the
[`digamma()`](https://rdrr.io/r/base/Special.html) and
[`trigamma()`](https://rdrr.io/r/base/Special.html) functions used in
the score and information steps. Also, the univariate optimization is
performed on `log(theta)`, while the Newton-Raphson approach is
performed on the original scale of `theta`.

## See also

[`MASS::theta.ml()`](https://rdrr.io/pkg/MASS/man/theta.md.html)

## Examples

``` r
library(MASS)

quine.nb <- glm.nb(Days ~ .^2, data = quine)
# estimate from MASS package
theta.ml(quine$Days, fitted(quine.nb), limit=200, eps=1e-6)
#> [1] 1.603641
#> attr(,"SE")
#> [1] 0.2138379

# faster for large-scale data
nb_theta(quine$Days, fitted(quine.nb))
#> [1] 1.603654
```
