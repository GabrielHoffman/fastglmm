# Variance Partitioning Analysis

Compute fraction of variance attributable to each variable in regression
model. Also interpretable as the intra-class correlation after
correcting for all other variables in the model.

## Usage

``` r
varpart(
  fit,
  ...,
  distr.method = c("trigamma", "lognormal"),
  lambda.method = c("parametric", "mean")
)

# S4 method for class 'fastlmm'
varpart(
  fit,
  ...,
  distr.method = c("trigamma", "lognormal"),
  lambda.method = c("parametric", "mean")
)

# S4 method for class 'fastglmm'
varpart(
  fit,
  ...,
  distr.method = c("trigamma", "lognormal"),
  lambda.method = c("parametric", "mean")
)

# S4 method for class 'glm'
varpart(
  fit,
  ...,
  distr.method = c("trigamma", "lognormal"),
  lambda.method = c("parametric", "mean")
)

# S4 method for class 'negbin'
varpart(
  fit,
  ...,
  distr.method = c("trigamma", "lognormal"),
  lambda.method = c("parametric", "mean")
)

# S4 method for class 'lm'
varpart(
  fit,
  ...,
  distr.method = c("trigamma", "lognormal"),
  lambda.method = c("parametric", "mean")
)

# S4 method for class 'merMod'
varpart(
  fit,
  ...,
  distr.method = c("trigamma", "lognormal"),
  lambda.method = c("parametric", "mean")
)
```

## Arguments

- fit:

  model fit

- ...:

  other arguments,

- distr.method:

  use either the `"lognormal"` or `"trigamma"` formulas from Nakagawa,
  et al. (2017)

- lambda.method:

  use either `"parametric"` or `"mean"` method to estimate the mean rate
  for count models

## Details

The coefficient of determination (i.e. R^2) is 1 - \[Residuals
fraction\]. This matches
[`performance::r2_nakagawa()`](https://easystats.github.io/performance/reference/r2_nakagawa.html)
and
[`performance::r2_mckelvey()`](https://easystats.github.io/performance/reference/r2_mckelvey.html),
except these use the `"lognormal"` method.

For count models, the distributional variance is a function of the mean
count rate. Following Eqn 5.8 of Nakagawa, et al. 2017, this can be
estimated using parameters of a model including only intercept and
random effect terms. But this requires refitting the model dropping the
rest of the fixed effects. Instead, computing the mean of the observed
counts is a fast approximation.

## References

Nakagawa, Johnson, Schielzeth. 2017. The coefficient of determination R2
and intra-class correlation coefficient from generalized linear
mixed-effects models revisited and expanded. J. R. Soc. Interface 14:
20170213.
[doi:10.1098/rsif.2017.0213](https://doi.org/10.1098/rsif.2017.0213)

Nakagawa, and Schielzeth. "A general and simple method for obtaining R2
from generalized linear mixed‐effects models." Methods in ecology and
evolution 4, no. 2 (2013): 133-142.
[doi:10.1111/j.2041-210x.2012.00261.x](https://doi.org/10.1111/j.2041-210x.2012.00261.x)

## Examples

``` r
library(MASS)
library(lme4)

fit = fastglmm(y ~ trt + I(week > 2) + (1 | ID),
  family = binomial(), data = bacteria)

varpart(fit)
#>         trt I(week > 2)          ID   Residuals 
#>  0.03695010  0.08324405  0.26113134  0.61867451 
```
