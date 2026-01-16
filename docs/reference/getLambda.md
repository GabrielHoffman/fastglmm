# Estimate baseline rate from count model

Estimate baseline rate from count model

Estimate baseline rate from count model

## Usage

``` r
getLambda(fit, method = c("parametric", "mean"), ...)

# S4 method for class 'fastlmm'
getLambda(fit, method = c("parametric", "mean"), ...)

# S4 method for class 'glm'
getLambda(fit, method = c("parametric", "mean"), ...)

# S4 method for class 'negbin'
getLambda(fit, method = c("parametric", "mean"), ...)

getLambdaFromNull(fit_null, ...)
```

## Arguments

- fit:

  model fit

- method:

  use either `"parametric"` or `"mean"` method to estimate the mean rate
  for count models

- ...:

  other args

- fit_null:

  model fit of null
