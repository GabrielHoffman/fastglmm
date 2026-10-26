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

## Value

model fit

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

if (FALSE) { # \dontrun{
# NB GLMM via Laplace approximation
fit1 <- glmer.nb(y ~ f1*f2 + (1|g), data=dd)
coef(summary(fit1))

# NB GLMM via PQL
fit2 <- fastglmm.nb(y ~ f1*f2 + (1|g), data=dd)
coef(summary(fit2))
} # }
```
