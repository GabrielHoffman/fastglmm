# Construct Design Matrices

Construct Design Matrices

## Usage

``` r
# S3 method for class 'fastlmm'
model.matrix(object, ...)
```

## Arguments

- object:

  regression model object

- ...:

  other args

## Examples

``` r
library(MASS)

# GLMM via PQL
fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

design <- model.matrix(fit)
head(design)
#>   (Intercept) trtdrug trtdrug+ I(week > 2)TRUE
#> 1           1       0        0               0
#> 2           1       0        0               0
#> 3           1       0        0               1
#> 4           1       0        0               1
#> 5           1       0        1               0
#> 6           1       0        1               0
```
