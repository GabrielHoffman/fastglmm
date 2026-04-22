# Plot Saturation Curve for Meff

Saturation curve for the effective number of independent measurements
per donor

## Usage

``` r
plotSaturation(fit, method = c("counts", "measurements"))
```

## Arguments

- fit:

  model fit from
  [`fastlmm()`](http://gabrielhoffman.github.io/fastglmm/reference/fastlmm.md)
  or
  [`fastglmm()`](http://gabrielhoffman.github.io/fastglmm/reference/fastglmm.md)

- method:

  use asymptotic values for `"counts"` or `"measurements"`

## Details

See
[`powerNBMM()`](http://gabrielhoffman.github.io/fastglmm/reference/powerNBMM.md)

## See also

[`meff()`](http://gabrielhoffman.github.io/fastglmm/reference/meff.md)
[`powerNBMM()`](http://gabrielhoffman.github.io/fastglmm/reference/powerNBMM.md)

## Examples

``` r
data(PsychAD)

# regression formula
form <- PTPRG ~ (1|SubID) + offset(log(libSize))

# fit NBMM in on PTPRG expression
fit <- fastglmm.nb(form, PsychAD)

# Saturation curve as mu increases
plotSaturation(fit, method="counts")


# Saturation curve as m increases
plotSaturation(fit, method="measurements")
```
