# Convert GLM family to string

Convert GLM family to string

## Usage

``` r
getFamilyString(family)
```

## Arguments

- family:

  description of the error distribution and link function to be used in
  the model, as a string or a function.

## Value

string value

## Examples

``` r
getFamilyString(gaussian)
#> [1] "gaussian/identity"

getFamilyString(gaussian())
#> [1] "gaussian/identity"

getFamilyString(poisson)
#> [1] "poisson/log"

getFamilyString(quasipoisson)
#> [1] "quasipoisson/log"

getFamilyString(quasibinomial)
#> [1] "quasibinomial/logit"

getFamilyString(quasibinomial("probit"))
#> [1] "quasibinomial/probit"

getFamilyString(binomial())
#> [1] "binomial/logit"

getFamilyString(binomial("probit"))
#> [1] "binomial/probit"

getFamilyString(MASS::negative.binomial(3))
#> [1] "nb:3"

getFamilyString(MASS::negative.binomial(NA))
#> [1] "nb"

getFamilyString("nb")
#> [1] "nb"
```
