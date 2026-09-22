# Is object a negative binomial model

Is object a negative binomial model

## Usage

``` r
isNB(x)
```

## Arguments

- x:

  family or model fit

## Value

TRUE for NB models

## Examples

``` r
data(PsychAD)

# regression formula
form <- PTPRG ~ (1|SubID) + offset(log(libSize))

# NB GLMM on PTPRG expression via PQL
fit <- fastglmm.nb(form, PsychAD)

isNB(fit)
#> [1] TRUE
```
