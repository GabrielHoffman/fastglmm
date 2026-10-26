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

  number of measurements per subject

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

## Value

object of class powerNBMM

## See also

[`meff()`](http://gabrielhoffman.github.io/fastglmm/reference/meff.md)

## Examples

``` r
data(PsychAD)

# regression formula
form <- PTPRG ~ (1|SubID) + offset(log(libSize))

# fit NBMM on PTPRG expression
# subset to reduce runtime
fit <- fastglmm.nb(form, PsychAD[seq(5000),])

# Power analysis
powerNBMM(fit = fit, beta = .1, sd_x= 0.5)
#>          N        f   lambda     power lambda.asymp.m lambda.asymp.mu
#> SubID 5000 1.311111 16.38888 0.9816171       15.72098        16.31384
#>       power.asymp.m power.asymp.mu  N.1   m     mu beta   sigSq.a  theta sd_x
#> SubID     0.9775187      0.9811948 5000 250 1.1454  0.1 0.7421914 0.2349  0.5
#>       alpha   kappa.m  kappa.mu n_measurements m.mean        rho    m.eff
#> SubID  0.05 0.9592463 0.9954213             20    250 0.08604907 11.14767
#>          m.max  fraction
#> SubID 11.62128 0.9592463
```
