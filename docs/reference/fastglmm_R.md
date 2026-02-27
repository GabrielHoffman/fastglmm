# Fit generalized linear mixed model via PQL

Fit generalized linear mixed model (GLMM) with a single random effect
using penalized quasi-likelihood (PQL)

## Usage

``` r
fastglmm_R(
  formula,
  data,
  family = gaussian(),
  weights = NULL,
  delta = NULL,
  delta.range = c(-10, 10),
  maxit = 100,
  tol = 1e-05,
  tol.eta = 1e-07,
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
fit = fastglmm_R(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)
coef(summary(fit))
#>                   Estimate Std. Error    df   t value     Pr(>|t|)
#> (Intercept)      3.4122672  0.5139248 101.1  6.639623 1.597027e-09
#> trtdrug         -1.2474256  0.6383952  46.9 -1.954002 5.668420e-02
#> trtdrug+        -0.7544013  0.6397064  49.9 -1.179293 2.438782e-01
#> I(week > 2)TRUE -1.6073742  0.3551299 184.3 -4.526159 1.075064e-05
```
