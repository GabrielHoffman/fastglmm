# Get link function given family

Get link function given family

## Usage

``` r
getLink(family)
```

## Arguments

- family:

  family function

## Value

link function for model family

## Examples

``` r
data(PsychAD)

# regression formula
form <- PTPRG ~ (1|SubID) + offset(log(libSize))

# NB GLMM on PTPRG expression via PQL
fit <- fastglmm.nb(form, PsychAD)

getLink(family(fit))
#> function (x, base = exp(1))  .Primitive("log")
```
