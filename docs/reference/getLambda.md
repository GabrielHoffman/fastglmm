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

## Value

scalar lambda value

scalar lambda value from the null model fit

## Details

For count models, the distributional variance is a function of the mean
count rate. Following Eqn 5.8 of Nakagawa, et al. 2017, this can be
estimated using parameters of a model including only intercept and
random effect terms. But this requires refitting the model dropping the
rest of the fixed effects. Instead, computing the mean of the observed
counts is a fast approximation.

## Examples

``` r
data(PsychAD)

# regression formula
form <- PTPRG ~ (1|SubID) + offset(log(libSize))

# NB GLMM on PTPRG expression via PQL
fit <- fastglmm.nb(form, PsychAD)

getLambda(fit)
#> Error: unable to find an inherited method for function ‘getLambda’ for signature ‘object = "fastglmm"’
data(PsychAD)

# regression formula
form <- PTPRG ~ (1|SubID) + offset(log(libSize))

# NB GLMM on PTPRG expression via PQL
fit <- fastglmm.nb(form, PsychAD)

getLambdaFromNull(fit)
#> Error: unable to find an inherited method for function ‘getLambdaFromNull’ for signature ‘fit_null = "fastglmm"’
```
