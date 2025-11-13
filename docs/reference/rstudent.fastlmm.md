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

## Examples

``` r
library(MASS)

# GLMM via PQL
fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

rstudent(fit)[1:3]
#> [1] 0.08309438 0.08309438 0.18953606
```
