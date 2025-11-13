# Simulate Responses

Simulate responses from model fit

## Usage

``` r
# S3 method for class 'fastlmm'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- object:

  model fit

- nsim:

  number of examples to simulate

- seed:

  random seed

- ...:

  other args, not used

## Examples

``` r
library(MASS)

# GLMM via PQL
fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

y.sim <- simulate(fit, 2)
head(y.sim)
#>      1_1 1_2
#> [1,]   1   1
#> [2,]   1   1
#> [3,]   1   1
#> [4,]   1   1
#> [5,]   1   1
#> [6,]   1   1
```
