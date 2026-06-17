# Power for Negative Binomial Mixed Model

Evaluate effective sample size and power for negative binomial mixed
model given parameter values

## Usage

``` r
powerNBMM(N, m, mu, sigSq.a, theta, beta = NA, sd_x = NA, alpha = 0.05, fit)
```

## Arguments

- N:

  number of subjects

- m:

  number of measurements pwer subject

- mu:

  mean read count

- sigSq.a:

  variance of random effect

- theta:

  negative binomial

- beta:

  effect size

- sd_x:

  standard deviation of target variable

- alpha:

  target false positive rate

- fit:

  model fit with
  [`fastlmm()`](http://gabrielhoffman.github.io/fastglmm/reference/fastlmm.md)
  or
  [`fastglmm()`](http://gabrielhoffman.github.io/fastglmm/reference/fastglmm.md)

## See also

[`meff()`](http://gabrielhoffman.github.io/fastglmm/reference/meff.md)

## Examples

``` r
data(PsychAD)

# regression formula
form <- PTPRG ~ (1|SubID) + offset(log(libSize))

# fit NBMM on PTPRG expression
fit <- fastglmm.nb(form, PsychAD)

# Power analysis
powerNBMM(fit = fit, beta = .1, sd_x= 0.5)
#>           N        f   lambda power lambda.asymp.m lambda.asymp.mu
#> SubID 60481 1.234716 186.6921     1        180.552        185.7482
#>       power.asymp.m power.asymp.mu   N.1        m       mu beta   sigSq.a
#> SubID             1              1 60481 202.9564 1.203237  0.1 0.7856561
#>        theta sd_x alpha   kappa.m  kappa.mu n_measurements   m.mean      rho
#> SubID 0.2445  0.5  0.05 0.9671114 0.9949439            298 202.9564 0.126551
#>          m.eff    m.max  fraction
#> SubID 7.642067 7.901952 0.9671114
```
