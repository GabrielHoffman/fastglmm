# Extract the Number of Observations from a Fit

Extract the Number of Observations from a Fit

## Usage

``` r
# S3 method for class 'fastlmm'
nobs(object, ...)
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

nobs(fit)
#> [1] 220
```
