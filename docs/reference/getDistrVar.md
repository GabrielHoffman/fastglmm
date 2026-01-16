# Distributional Variance

Compute distributional variance from the model fit

## Usage

``` r
getDistrVar(
  fit,
  fit_null,
  method = c("trigamma", "lognormal"),
  lambda.method = c("parametric", "mean")
)
```

## Arguments

- fit:

  model fit

- method:

  use either the `"lognormal"` or `"trigamma"` formulas from Nakagawa,
  et al. (2017)

- lambda.method:

  use either `"parametric"` or `"mean"` method to estimate the mean rate
  for count models

## Details

In generalized linear (mixed) models, the link function contributes to
the coefficient of determination (Nakagawa, et al., 2012, 2017; McKelvey
and Zavoina, 1975).

1 - Residuals gives the R2 values from
`performance::r2_nakagawa(..., approximation="trigamma")`. Using
[`performance::r2_mckelvey()`](https://easystats.github.io/performance/reference/r2_mckelvey.html)
use the "lognormal" approximation

For count models, the distributional variance is a function of the mean
count rate. Following Eqn 5.8 of Nakagawa, et al. 2017, this can be
estimated using parameters of a model including only intercept and
random effect terms. But this requires refitting the model dropping the
rest of the fixed effects. Instead, computing the mean of the observed
counts is a fast approximation.

## References

Nakagawa, Johnson, Schielzeth. 2017. The coefficient of determination R2
and intra-class correlation coefficient from generalized linear
mixed-effects models revisited and expanded. J. R. Soc. Interface 14:
20170213.
[doi:10.1098/rsif.2017.0213](https://doi.org/10.1098/rsif.2017.0213)

Nakagawa, Shinichi, and Holger Schielzeth. "A general and simple method
for obtaining R2 from generalized linear mixed‐effects models." Methods
in ecology and evolution 4, no. 2 (2013): 133-142.
[doi:10.1111/j.2041-210x.2012.00261.x](https://doi.org/10.1111/j.2041-210x.2012.00261.x)

McKelvey, R., Zavoina, W. (1975), "A Statistical Model for the Analysis
of Ordinal Level Dependent Variables", Journal of Mathematical Sociology
4, S.103–120.
