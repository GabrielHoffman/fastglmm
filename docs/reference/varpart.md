# Variance Partitioning Analysis

Compute fraction of variance attributable to each variable in regression
model. Also interpretable as the intra-class correlation after
correcting for all other variables in the model.

## Usage

``` r
varpart(fit, ..., distr.method = c("trigamma", "lognormal"))

# S4 method for class 'fastlmm'
varpart(fit, ..., distr.method = c("trigamma", "lognormal"))

# S4 method for class 'fastglmm'
varpart(fit, ..., distr.method = c("trigamma", "lognormal"))

# S4 method for class 'glm'
varpart(fit, ..., distr.method = c("trigamma", "lognormal"))

# S4 method for class 'lm'
varpart(fit, ..., distr.method = c("trigamma", "lognormal"))

# S4 method for class 'merMod'
varpart(fit, ..., distr.method = c("trigamma", "lognormal"))
```

## Arguments

- fit:

  model fit

- ...:

  other arguments, not used here

- distr.method:

  use either the `"lognormal"` or `"trigamma"` formulas from Nakagawa,
  et al. (2017)

## Details

The coefficient of determination (i.e. R^2) is 1 - \[Residuals
fraction\]. This matches `performance::r2_nakagawa()` and
`performance::r2_mckelvey()`, except these use the `"lognormal"` method.

## References

Nakagawa, Johnson, Schielzeth. 2017. The coefficient of determination R2
and intra-class correlation coefficient from generalized linear
mixed-effects models revisited and expanded. J. R. Soc. Interface 14:
20170213.
[doi:10.1098/rsif.2017.0213](https://doi.org/10.1098/rsif.2017.0213)

Nakagawa, Shinichi, and Holger Schielzeth. "A general and simple method
for obtaining R2 from generalized linear mixed‐effects models." Methods
in ecology and evolution 4, no. 2 (2013): 133-142.
[doi:10.1111/j.2041-210x.2012.00261.x](https://doi.org/10.1111/j.2041-210x.2012.00261.x)

## Examples

``` r
library(MASS)
library(lme4)

fit = fastglmm(y ~ trt + I(week > 2) + (1 | ID),
  family = binomial(), data = bacteria)

varpart(fit)
#>         trt I(week > 2)          ID   Residuals 
#>  0.03696502  0.08327439  0.26104630  0.61871429 
```
