# Extract Model Fitted Values

Extract Model Fitted Values

## Usage

``` r
# S3 method for class 'fastlmm'
fitted(object, ...)

# S3 method for class 'fastglmm'
fitted(object, ...)
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

fitted(fit)[1:3]
#> [1] 0.9828584 0.9828584 0.9199276
```
