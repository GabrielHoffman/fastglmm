# Residual Degrees-of-Freedom

Residual degrees-of-freedom is the trace of the residual hat matrix

## Usage

``` r
# S3 method for class 'fastlmm'
df.residual(object, ...)
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

df.residual(fit)
#> [1] 189.1795
```
