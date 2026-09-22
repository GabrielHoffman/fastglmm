# Effective Degrees-of-Freedom of Model Fit

Effective degrees-of-freedom of model fit is the trace of the hat matrix
mapping from observed to predicted response

## Usage

``` r
edf(object, ...)
```

## Arguments

- object:

  fitted model

- ...:

  other args, not used

## Value

effective degrees of freedom

## Examples

``` r
library(MASS)

# GLMM via PQL
fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

edf(fit)
#> [1] 30.82029
```
