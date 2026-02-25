# Calculate Variance-Covariance Matrix for a Fitted Model Object

Calculate Variance-Covariance Matrix for a Fitted Model Object

## Usage

``` r
# S3 method for class 'fastlmm'
vcov(object, ...)
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

vcov(fit)
#>                 (Intercept)     trtdrug      trtdrug+ I(week > 2)TRUE
#> (Intercept)       0.4294063 -0.31905938 -0.3052379240   -0.1592879082
#> trtdrug          -0.3190594  0.66270396  0.3055753266    0.0172978631
#> trtdrug+         -0.3052379  0.30557533  0.6654168084   -0.0004870472
#> I(week > 2)TRUE  -0.1592879  0.01729786 -0.0004870472    0.2049654054
```
