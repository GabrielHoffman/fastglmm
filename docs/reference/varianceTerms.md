# Array of variance component estimates

Array of variance component estimates

## Usage

``` r
varianceTerms(object, ...)

# S4 method for class 'lm'
varianceTerms(object, ...)

# S4 method for class 'glm'
varianceTerms(object, ...)

# S4 method for class 'negbin'
varianceTerms(object, ...)

# S4 method for class 'merMod'
varianceTerms(object, ...)

# S4 method for class 'glmmTMB'
varianceTerms(object, ...)

# S4 method for class 'fastlmm'
varianceTerms(object, ...)
```

## Arguments

- object:

  model fit

- ...:

  other args

## Value

array of variance component estimates

## Examples

``` r
data(PsychAD)

# regression formula
form <- PTPRG ~ (1|SubID) + offset(log(libSize))

# NB GLMM on PTPRG expression via PQL
fit <- fastglmm.nb(form, PsychAD)

varianceTerms(fit)
#>     SubID 
#> 0.7856561 
```
