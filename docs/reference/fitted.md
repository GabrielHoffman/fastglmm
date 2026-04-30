# Extract Model Fitted Values

Extract Model Fitted Values

## Usage

``` r
# S3 method for class 'fastlmm'
fitted(object, ..., newdata = NULL)

# S3 method for class 'fastglmm'
fitted(object, ..., newdata = NULL)
```

## Arguments

- object:

  fitted model of class `fastlmm`

- ...:

  other args, not used

- newdata:

  optionally, a data.frame in which to look for variables with which to
  predict. If omitted, the fitted predictors are used.

## Examples

``` r
library(MASS)

# GLMM via PQL
fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

fitted(fit)[1:3]
#> [1] 0.9828584 0.9828584 0.9199276
```
