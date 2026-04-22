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

- nthreads:

  number of threads

## Examples

``` r
library(MASS)
data(PsychAD)

# regression formula
form <- PTPRG ~ (1|SubID) + offset(log(libSize))

# NB GLMM on PTPRG expression via PQL
fit1 <- fastglmm.nb(form, data)
#> Error in as.data.frame.default(data, optional = TRUE): cannot coerce class ‘"function"’ to a data.frame
coef(summary(fit1))
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'summary': object 'fit1' not found

# NB GLMM via Laplace approximation
fit2 <- lme4::glmer.nb(form, data)
#> Error in list2env(data): first argument must be a named list
coef(summary(fit2))
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'summary': object 'fit2' not found
```
