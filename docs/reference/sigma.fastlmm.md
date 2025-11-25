# Extract Residual Standard Deviation

Extract Residual Standard Deviation

## Usage

``` r
# S3 method for class 'fastlmm'
sigma(object, ...)
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

sigma(fit)
#> [1] 2.262139
```
