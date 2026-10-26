# Fit Generalized Linear Mixed Model via PQL

Fit generalized linear mixed model (GLMM) with a single random effect
using penalized quasi-likelihood (PQL)

## Usage

``` r
fastglmm(
  formula,
  data,
  family = gaussian(),
  weights = NULL,
  delta = NULL,
  delta.range = c(-10, 10),
  maxit = 1000,
  tol = 0.001,
  tol.eta = 0.001,
  doCoxReid = TRUE,
  lambda = 0,
  nthreads = 6,
  REML = TRUE
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

- family:

  a description of the error distribution and link function to be used
  in the model.

- weights:

  an optional vector of prior weights with a value for each sample.

- delta:

  if `NULL` estimate delta, if value is given uses this fixed value

- delta.range:

  min and max values (in log space), of the search space for delta to
  fit the random effect

- maxit:

  max number of PQL iterations

- tol:

  convergence criterion for the 1D search of the delta space

- tol.eta:

  convergence criterion `eta` in the PQL iteration

- doCoxReid:

  use Cox-Reid correction for estimating theta in negative binomial
  model

- lambda:

  ridge shrinkage parameter

- nthreads:

  number of threads

- REML:

  use REML or ML estimates (default to `FALSE`)

## Value

model fit object of class fastglmm

## Details

Parameters are estimated using penalized quasi-likelihood (PQL) (Breslow
and Clayton, 1993). See Maestrini, et al (2026) for an examination of
REML in the context of generalized linear mixed models

## References

Breslow, N.E. and Clayton, D.G. (1993). Approximate inference in
generalized linear mixed models. Journal of the American statistical
Association, 88(421), pp.9-25.
[doi:10.1080/01621459.1993.10594284](https://doi.org/10.1080/01621459.1993.10594284)

Maestrini, L., Hui, F.K. and Welsh, A.H., (2026). Restricted maximum
likelihood estimation in generalized linear mixed models. Statistical
Science, 41(3) 615-636
[doi:10.1214/25-STS995](https://doi.org/10.1214/25-STS995)

## Examples

``` r
library(MASS)
library(lme4)

# GLMM via Laplace approximation
fit = glmer(y ~ trt + I(week > 2) + (1 | ID),
  family = binomial(), data = bacteria)
coef(summary(fit))
#>                   Estimate Std. Error   z value     Pr(>|z|)
#> (Intercept)      3.5479411  0.6957868  5.099178 3.411311e-07
#> trtdrug         -1.3666520  0.6769898 -2.018719 4.351647e-02
#> trtdrug+        -0.7826412  0.6831225 -1.145682 2.519267e-01
#> I(week > 2)TRUE -1.5984897  0.4759391 -3.358601 7.833796e-04

# GLMM via PQL
fit = glmmPQL(y ~ trt + I(week > 2), random = ~ 1 | ID,
  family = binomial, data = bacteria, verbose = FALSE)
coef(summary(fit))
#>                      Value Std.Error  DF   t-value      p-value
#> (Intercept)      3.4120140 0.5185033 169  6.580506 5.646616e-10
#> trtdrug         -1.2473553 0.6440635  47 -1.936696 5.880736e-02
#> trtdrug+        -0.7543273 0.6453978  47 -1.168779 2.483867e-01
#> I(week > 2)TRUE -1.6072570 0.3583379 169 -4.485311 1.340324e-05

# GLMM via PQL
fit = fastglmm(y ~ trt + I(week > 2) + (1 | ID),
  family = binomial(), data = bacteria)
coef(summary(fit))
#>                  Estimate Std. Error     df    t value     Pr(>|t|)
#> (Intercept)      3.475641  0.6843895 223.72  5.0784541 8.010499e-07
#> trtdrug         -1.265157  0.8610477 105.31 -1.4693228 1.447264e-01
#> trtdrug+        -0.773049  0.8612327 111.54 -0.8976075 3.713283e-01
#> I(week > 2)TRUE -1.636336  0.4599847 490.21 -3.5573703 4.109830e-04
```
