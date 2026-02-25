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

## Examples

``` r
library(MASS)

# GLMM via PQL
fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

data <- model.frame(fit)
head(data)
#>   y     trt I(week > 2)  ID (weights)
#> 1 y placebo       FALSE X01 0.1415670
#> 2 y placebo       FALSE X01 0.1415670
#> 3 y placebo        TRUE X01 0.6189372
#> 4 y placebo        TRUE X01 0.6189372
#> 5 y   drug+       FALSE X02 0.7217012
#> 6 y   drug+       FALSE X02 0.7217012
```
