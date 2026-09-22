# Get theta from NB model

Get theta from NB model

## Usage

``` r
getTheta(fit)
```

## Arguments

- fit:

  model fit

## Value

theta for NB models, else NA

## Examples

``` r
data(PsychAD)

# regression formula
form <- PTPRG ~ (1|SubID) + offset(log(libSize))

# NB GLMM on PTPRG expression via PQL
fit <- fastglmm.nb(form, PsychAD)

getTheta(fit)
#> Error: unable to find an inherited method for function ‘getTheta’ for signature ‘object = "fastglmm"’
```
