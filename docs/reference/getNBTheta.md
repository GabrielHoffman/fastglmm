# Get theta from NB model

Get theta from NB model

## Usage

``` r
getNBTheta(fit)
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

getNBTheta(fit)
#> [1] 0.2445
```
