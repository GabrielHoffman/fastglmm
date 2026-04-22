# Effective Number of Independent Measurements per Subject

Compute the effective number of independent measurements per subject for
a generalized linear mixed model

## Usage

``` r
meff(fit)
```

## Arguments

- fit:

  model fit from
  [`fastlmm()`](http://gabrielhoffman.github.io/fastglmm/reference/fastlmm.md)
  or
  [`fastglmm()`](http://gabrielhoffman.github.io/fastglmm/reference/fastglmm.md)

## Value

- `m.mean`: mean number of measurements per subject

- `rho`: intra-class correlation within subjects

- `m.eff`: effective number of independent measurements per subject

- `m.max`: maximum value of m.eff for m -\> Inf

- `fraction`: m.eff / m.max

## Details

In a repeated measures model with multiple correlated measurements per
subject, Lui and Liang (1997) define a formula for the effect number of
\_independent\_ measurements per subject (m.eff): \$\$ m\_\text{eff} = m
/ (1+(m-1)\rho) \$\$ where \\m\\ is the mean number of measurements per
subjectd, and \\\rho\\ is the intra-class correlation indicating the
correlation between measurements within the same subject.

Increasing \\m\\ has diminishing returns and as \\m\\ increases,
\\m\_\text{eff}\\ converges to \\1/\rho\\.

## References

Liu, G. and Liang, K.Y., 1997. Sample size calculations for studies with
correlated observations. Biometrics, pp.937-947.

## See also

`plotMeff()`

## Examples

``` r
data(PsychAD)

# regression formula
form <- PTPRG ~ (1|SubID) + offset(log(libSize))

# fit NBMM in on PTPRG expression
fit <- fastglmm.nb(form, PsychAD)

# effective number of independent measurements per subject
meff(fit)
#>     m.mean        rho    m.eff    m.max  fraction
#> 1 202.9564 0.03465127 25.37575 28.85897 0.8793019
```
