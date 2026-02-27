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
#> (Intercept)  1.6500844  0.1298536 12.707264 5.388716e-37
#> f12          0.7671505  0.1815543  4.225460 2.384528e-05
#> f13          1.0114682  0.1811596  5.583300 2.359976e-08
#> f2B          1.5124145  0.1805961  8.374570 5.542585e-17
#> f12:f2B     -0.6150597  0.2538182 -2.423229 1.538322e-02
#> f13:f2B     -0.6103981  0.2534190 -2.408652 1.601155e-02

# NB GLMM via PQL
fit2 <- fastglmm.nb(y ~ f1*f2 + (1|g), data=dd)
coef(summary(fit2))
#>               Estimate Std. Error    df   t value     Pr(>|t|)
#> (Intercept)  1.6500365  0.1298459 272.1 12.707655 2.322370e-29
#> f12          0.7671854  0.1814985 802.6  4.226952 2.640988e-05
#> f13          1.0115619  0.1811065 802.6  5.585454 3.191837e-08
#> f2B          1.5122439  0.1805473 802.6  8.375890 2.433819e-16
#> f12:f2B     -0.6149961  0.2537603 802.6 -2.423531 1.559056e-02
#> f13:f2B     -0.6103809  0.2533639 802.6 -2.409108 1.621583e-02
```
