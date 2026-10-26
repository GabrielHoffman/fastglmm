# Extract Model Coefficients

Extract Model Coefficients

## Usage

``` r
# S3 method for class 'fastlmm'
coef(object, ...)
```

## Arguments

- object:

  fitted model of class `fastlmm`

- ...:

  other args, not used

## Value

coefficient values

## Examples

``` r
library(MASS)

# GLMM via PQL
fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

coef(fit)
#>     (Intercept)         trtdrug        trtdrug+ I(week > 2)TRUE 
#>        3.475641       -1.265157       -0.773049       -1.636336 
```
