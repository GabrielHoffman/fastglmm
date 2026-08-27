# Fit Negative Binomial Mixed Model via PQL

Fit negative binomial mixed model (GLMM) with a single random effect
using penalized quasi-likelihood (PQL)

## Usage

``` r
fastglmm.nb(
  formula,
  data,
  weights = NULL,
  maxit = 100,
  tol = 0.001,
  tol.eta = 0.001,
  doCoxReid = nrow(data) < 1000,
  lambda = 0,
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

- doCoxReid:

  use Cox-Reid correction for estimating theta in negative binomial
  model

- lambda:

  ridge shrinkage parameter

- nthreads:

  number of threads

## Examples

``` r
library(MASS)
data(PsychAD)

# regression formula
form <- PTPRG ~ (1|SubID) + offset(log(libSize))

# NB GLMM on PTPRG expression via PQL
fit1 <- fastglmm.nb(form, PsychAD)
coef(summary(fit1))
#>              Estimate Std. Error    df   t value      Pr(>|t|)
#> (Intercept) -8.638731 0.04585474 148.2 -188.3934 2.647753e-178

# NB GLMM via Laplace approximation
# fit2 <- lme4::glmer.nb(form, PsychAD)
# coef(summary(fit2))
```
