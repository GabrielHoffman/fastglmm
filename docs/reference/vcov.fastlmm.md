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
#>                 (Intercept)    trtdrug      trtdrug+ I(week > 2)TRUE
#> (Intercept)      0.26423415 -0.1963324 -0.1878274432   -0.0980174420
#> trtdrug         -0.19633245  0.4077933  0.1880350631    0.0106441996
#> trtdrug+        -0.18782744  0.1880351  0.4094626779   -0.0002997034
#> I(week > 2)TRUE -0.09801744  0.0106442 -0.0002997034    0.1261249832
```
