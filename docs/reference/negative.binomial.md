# Family function for Negative Binomial GLMs

Specifies the information required to fit a Negative Binomial
generalized linear model, with known `theta` parameter.

## Usage

``` r
negative.binomial(theta = stop("'theta' must be specified"), link = "log")
```

## Arguments

- theta:

  The known value of the additional parameter, `theta`. If set to `NA`,
  [`fastglmm()`](http://gabrielhoffman.github.io/fastglmm/reference/fastglmm-pkg.md)
  will estimate the parameter from the data

- link:

  The link function, as a character string, name or one-element
  character vector specifying one of `log`, `sqrt` or `identity`, or an
  object of class `"link-glm"`.
