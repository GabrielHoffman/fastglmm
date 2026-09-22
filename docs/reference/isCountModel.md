# Is object a count model

Is object a count model

## Usage

``` r
isCountModel(x)
```

## Arguments

- x:

  family or model fit

## Value

TRUE for poisson, quasipoisson or NB models

## Examples

``` r
data(PsychAD)

# regression formula
form <- PTPRG ~ (1|SubID) + offset(log(libSize))

# NB GLMM on PTPRG expression via PQL
fit <- fastglmm.nb(form, PsychAD)

isCountModel(fit)
#> [1] TRUE
```
