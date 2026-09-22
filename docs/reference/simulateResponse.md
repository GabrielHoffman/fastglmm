# Simulate Responses

Simulate responses from model fit

## Usage

``` r
simulateResponse(mu, nsim, family, sd, theta)
```

## Arguments

- mu:

  condition on this systematic component

- nsim:

  number of examples to simulate

- family:

  regression family

- sd:

  standard deviation from model fit

- theta:

  NB theta

## Value

matrix of responses simulated from the model

## Examples

``` r
data(PsychAD)

# regression formula
form <- PTPRG ~ (1|SubID) + offset(log(libSize))

# NB GLMM on PTPRG expression via PQL
fit <- fastglmm.nb(form, PsychAD)

Y <- simulateResponse(fitted(fit), 10, "nb", sigma(fit), getNBTheta(fit))

Y[1:2, 1:3]
#>      [,1] [,2] [,3]
#> [1,]    2    4    0
#> [2,]    0    0    0
```
