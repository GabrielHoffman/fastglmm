# Dispersion parameter phi for quasi-likelihood

Extract dispersion parameter phi for quasi-likelihood

## Usage

``` r
dispersion(object)

# S4 method for class 'fastlmm'
dispersion(object)

# S4 method for class 'glm'
dispersion(object)

# S4 method for class 'negbin'
dispersion(object)

# S4 method for class 'glmmTMB'
dispersion(object)

# S4 method for class 'glmerMod'
dispersion(object)
```

## Arguments

- object:

  model fit

## Examples

``` r
library(MASS)
data(PsychAD)

# regression formula
form <- PTPRG ~ (1|SubID) + offset(log(libSize))

# NB GLMM on PTPRG expression via PQL
fit <- fastglmm.nb(form, PsychAD)

dispersion(fit)
#> [1] 1
```
