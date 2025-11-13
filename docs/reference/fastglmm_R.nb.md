# Fit negative binomial mixed model via PQL

Fit negative binomial mixed model (GLMM) with a single random effect
using penalized quasi-likelihood (PQL)

## Usage

``` r
fastglmm_R.nb(
  formula,
  data,
  weights,
  maxit = 100,
  tol = .Machine$double.eps^0.5,
  tol.eta = .Machine$double.eps^0.5,
  init.fit = NULL,
  init = c("lm", "glm"),
  nthreads = 6
)
```

## Arguments

- formula:

  a two-sided linear formula object describing both the fixed-effects
  and random-effects part of the model, with the response on the left of
  a `~` operator and the terms, separated by `+` operators, on the
  right. Random-effects terms are distinguished by vertical bars (`|`)
  separating expressions for design matrices from grouping factors.

- data:

  an optional data frame containing the variables named in

- weights:

  an optional vector of prior weights with a value for each sample.

- maxit:

  max number of NB iterations

- tol:

  convergence criterion for the 1D search of the delta space

- tol.eta:

  convergence criterion `eta` in the PQL iteration

- init.fit:

  `fastglmm` object to initialize parameters

- init:

  `c("lm", "glm")` method to initialize `eta` values

- nthreads:

  number of threads

## Examples

``` r
library(MASS)
library(lme4)

set.seed(101)
dd <- expand.grid(f1 = factor(1:3),
               f2 = LETTERS[1:2], g=factor(1:9), rep=1:15,
       KEEP.OUT.ATTRS=FALSE)
mu <- 5*(-4 + with(dd, as.integer(f1) + 4*as.numeric(f2)))
dd$y <- rnbinom(nrow(dd), mu = mu, size = 0.5)

# NB GLMM via Laplace approximation
fit1 <- glmer.nb(y ~ f1*f2 + (1|g), data=dd)
#> boundary (singular) fit: see help('isSingular')
coef(summary(fit1))
#>               Estimate Std. Error   z value     Pr(>|z|)
#> (Intercept)  1.6500844  0.1298550 12.707129 5.398042e-37
#> f12          0.7671505  0.1815542  4.225462 2.384511e-05
#> f13          1.0114682  0.1811614  5.583242 2.360756e-08
#> f2B          1.5124145  0.1805975  8.374503 5.545753e-17
#> f12:f2B     -0.6150597  0.2538168 -2.423243 1.538264e-02
#> f13:f2B     -0.6103981  0.2534214 -2.408629 1.601258e-02

# NB GLMM via PQL
fit2 <- fastglmm.nb(y ~ f1*f2 + (1|g), data=dd)
coef(summary(fit2))
#>               Estimate Std. Error   z value     Pr(>|z|)
#> (Intercept)  1.6500367  0.1298045 12.711707 5.091062e-37
#> f12          0.7671853  0.1814267  4.228624 2.351245e-05
#> f13          1.0115617  0.1810321  5.587747 2.300346e-08
#> f2B          1.5122440  0.1804692  8.379511 5.314793e-17
#> f12:f2B     -0.6149960  0.2536396 -2.424684 1.532171e-02
#> f13:f2B     -0.6103808  0.2532406 -2.410281 1.594025e-02
```
