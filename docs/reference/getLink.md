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

getLink(fit)
#> Error in getFamilyString(family): family must be a string or function
```
