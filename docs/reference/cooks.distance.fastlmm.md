# Cook's Distance Metric

Cook's Distance Metric

## Usage

``` r
# S3 method for class 'fastlmm'
cooks.distance(model, ...)
```

## Arguments

- model:

  fitted model of class `fastlmm`

- ...:

  other args, not used

## Examples

``` r
library(MASS)

# GLMM via PQL
fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

cooks.distance(fit)[1:3]
#> [1] 3.589912e-05 3.589912e-05 9.738783e-04
```
