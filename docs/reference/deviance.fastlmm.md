# Model Deviance

Model Deviance

## Usage

``` r
# S3 method for class 'fastlmm'
deviance(object, ...)
```

## Arguments

- object:

  fitted model of class `fastlmm`

- ...:

  other args, not used

## Examples

``` r
library(MASS)

# GLMM via PQL
fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

deviance(fit)
#> [1] 1101.846
```
