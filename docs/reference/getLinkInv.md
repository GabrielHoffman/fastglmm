# Get inverse link function given family

Get inverse link function given family

## Usage

``` r
getLinkInv(family)
```

## Arguments

- family:

  family function

## Value

inverse link function for model family

## Examples

``` r
data(PsychAD)

# regression formula
form <- PTPRG ~ (1|SubID) + offset(log(libSize))

# NB GLMM on PTPRG expression via PQL
fit <- fastglmm.nb(form, PsychAD)

getLinkInv(family(fit))
#> function (x)  .Primitive("exp")
```
