# Extract AIC from a Fitted Model

Extract AIC from a Fitted Model

## Usage

``` r
# S3 method for class 'fastlmm'
extractAIC(fit, scale = 0, k = 2, ...)
```

## Arguments

- fit:

  fitted model of class `fastlmm`

- scale:

  not used

- k:

  numeric specifying the `weight` of the \_equivalent degrees of
  freedom\_ (=: `edf`) part in the AIC formula.

- ...:

  other args, not used

## Examples

``` r
library(MASS)

# GLMM via PQL
fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

extractAIC(fit)
#> [1]   30.82054 1163.45215
```
