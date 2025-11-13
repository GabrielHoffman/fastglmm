# Extract Log-Likelihood

Extract Log-Likelihood

## Usage

``` r
# S3 method for class 'fastlmm'
logLik(object, ...)
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

logLik(fit)
#> 'log Lik.' -550.9055 (df=30.82054)
```
