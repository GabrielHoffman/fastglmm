# Return Diagonals of Hat Matrix

Return diagonals of the hat matrix

## Usage

``` r
# S3 method for class 'fastlmm'
hatvalues(model, ...)
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

hatvalues(fit)[1:3]
#>          1          2          3 
#> 0.03886899 0.03886899 0.16102200 
```
