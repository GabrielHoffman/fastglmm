# Extracting the Model Frame from a Fit

Extracting the Model Frame from a Fit

## Usage

``` r
# S3 method for class 'fastlmm'
model.frame(formula, ...)
```

## Arguments

- formula:

  regression fit

- ...:

  other args

## Value

model frame

## Examples

``` r
library(MASS)

# GLMM via PQL
fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

data <- model.frame(fit)
head(data)
#>   y     trt I(week > 2)  ID (weights)
#> 1 y placebo       FALSE X01 0.1303500
#> 2 y placebo       FALSE X01 0.1303500
#> 3 y placebo        TRUE X01 0.5913423
#> 4 y placebo        TRUE X01 0.5913423
#> 5 y   drug+       FALSE X02 0.7284608
#> 6 y   drug+       FALSE X02 0.7284608
```
