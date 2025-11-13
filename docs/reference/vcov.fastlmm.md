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
#> (Intercept)      0.26412097 -0.19622922 -0.1877251075   -0.0980071676
#> trtdrug         -0.19622922  0.40755251  0.1879333561    0.0106426435
#> trtdrug+        -0.18772511  0.18793336  0.4092282213   -0.0003006235
#> I(week > 2)TRUE -0.09800717  0.01064264 -0.0003006235    0.1261176728
```
