# Simulate Responses

Simulate responses from model fit, conditional on the random effects

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

## See also

[`lme4::simulate.merMod()`](https://rdrr.io/pkg/lme4/man/simulate.merMod.html)

## Examples

``` r
library(MASS)

# GLMM via PQL
fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

y.sim <- simulate(fit, 2)
head(y.sim)
#>   sim_1 sim_2
#> 1     1     1
#> 2     1     1
#> 3     1     1
#> 4     1     1
#> 5     1     1
#> 6     1     1
```
