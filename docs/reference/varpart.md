# Variance Partitioning Analysis

Compute fraction of variance attributable to each variable in regression
model. Also interpretable as the intra-class correlation after
correcting for all other variables in the model.

## Usage

``` r
varpart(
  fit,
  method = c("exact", "approximate", "trigamma", "lognormal", "delta"),
  pseudocount = 1,
  p.tail = 1e-04,
  ...
)

# S4 method for class 'fastlmm'
varpart(
  fit,
  method = c("exact", "approximate", "trigamma", "lognormal", "delta"),
  pseudocount = 1,
  p.tail = 1e-04,
  ...
)

# S4 method for class 'fastglmm'
varpart(
  fit,
  method = c("exact", "approximate", "trigamma", "lognormal", "delta"),
  pseudocount = 1,
  p.tail = 1e-04,
  ...
)

# S4 method for class 'glm'
varpart(
  fit,
  method = c("exact", "approximate", "trigamma", "lognormal", "delta"),
  pseudocount = 1,
  p.tail = 1e-04,
  ...
)

# S4 method for class 'negbin'
varpart(
  fit,
  method = c("exact", "approximate", "trigamma", "lognormal", "delta"),
  pseudocount = 1,
  p.tail = 1e-04,
  ...
)

# S4 method for class 'lm'
varpart(
  fit,
  method = c("exact", "approximate", "trigamma", "lognormal", "delta"),
  pseudocount = 1,
  p.tail = 1e-04,
  ...
)

# S4 method for class 'merMod'
varpart(
  fit,
  method = c("exact", "approximate", "trigamma", "lognormal", "delta"),
  pseudocount = 1,
  p.tail = 1e-04,
  ...
)

# S4 method for class 'glmmTMB'
varpart(
  fit,
  method = c("exact", "approximate", "trigamma", "lognormal", "delta"),
  pseudocount = 1,
  p.tail = 1e-04,
  ...
)
```

## Arguments

- fit:

  regression model fit

- method:

  select method for count models: `"exact"` or `"approximate"` from the
  current work, or `"trigamma"`, `"lognormal"` or `"delta"` formulas
  from Nakagawa, et al. (2017)

- pseudocount:

  pseudocount used for `"exact"` and `"approximate"` methods for count
  models

- p.tail:

  probability threashold for evaluating expectations for `"exact"`
  methods for count models

- ...:

  other arguments, passed to `vpOther()` or `vpCounts()`

## Details

For linear model, variance fractions are computed based on the sum of
squares explained by each component. For the linear mixed model, the
variance fractions are computed by variance component estimates for
random effects and sum of squares for fixed effects.

For a generalized linear model, the variance fraction also includes the
contribution of the link function so that fractions are reported on the
linear (i.e. link) scale rather than the observed (i.e. response) scale.
For linear regression with an identity link, fractions are the same on
both scales. But for logit or probit links, the fractions are not well
defined on the observed scale due to the transformation imposed by the
link function.

The variance implied by the link function is the variance of the
corresponding distribution (Nakagawa, et al. 2013, 2017)

logit -\> logistic distribution -\> variance is \\\pi^\frac{2}{3}\\

probit -\> standard normal distribution -\> variance is 1

For count models, Nakagawa, et al. (2013, 2017) propose a large-count
approximation. Instead, we use an exact method described in Hoffman, et
al (2026).

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

Hoffman, et al. Partitioning gene expression variance using count
models. In prep.

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
