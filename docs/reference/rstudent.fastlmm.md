# Extract Studentized Residuals

Extract Studentized Residuals

## Usage

``` r
# S3 method for class 'fastlmm'
rstudent(model, ...)
```

## Arguments

- model:

  fitted model of class `fastlmm`

- ...:

  other args, not used

## Value

standardized residuals

## Examples

``` r
library(MASS)

# GLMM via PQL
fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

rstudent(fit)[1:3]
#> [1] 0.07810746 0.07810746 0.18136513
```
