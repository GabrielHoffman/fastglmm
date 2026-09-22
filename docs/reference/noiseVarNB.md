# Evaluate NB variance when lambda is Inf

Evaluate NB variance when lambda is Inf

## Usage

``` r
noiseVarNB(theta, method = c("trigamma", "lognormal", "delta"))
```

## Arguments

- theta:

  NB overdispersion parameter

- method:

  approximation method

## Value

scalar value the the NB variance

## Examples

``` r
data(PsychAD)

# regression formula
form <- PTPRG ~ (1|SubID) + offset(log(libSize))

# NB GLMM on PTPRG expression via PQL
fit <- fastglmm.nb(form, PsychAD)

noiseVarNB(getNBTheta(fit))
#> [1] 17.93261
```
