# Denominator Degrees of Freedom

Denominator degrees of freedom using Satterthwaite method

## Usage

``` r
.ddf(delta, hessian.vc, A.sat, B.sat, V, L)
```

## Arguments

- delta:

  ratio of variance components

- hessian.vc:

  hession

- A.sat:

  A matrix

- B.sat:

  B matrix

- V:

  variance-covariance matrix

- L:

  matrix of coefficient contrasts, one per \_row\_

## Value

array, denominator degrees of freedom for each contrast (i.e. \_row\_)

## Examples

``` r
library(MASS)

# GLMM via PQL
fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

# denominator degrees of freedom 
# used for hypothesis testing below
ddf(fit)
#> [1] 267.0785 123.8255 131.7683 485.7446
```
