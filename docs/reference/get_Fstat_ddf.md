# DDF for joint test

DDF for joint test by summarizing ddf for each coefficient

## Usage

``` r
get_Fstat_ddf(nu, tol = 1e-08)
```

## Arguments

- nu:

  array of ddf values

- tol:

  tolerance

## Value

scalar value

## See also

`lmerTest:::get_Fstat_ddf()`

## Examples

``` r
get_Fstat_ddf(1:10)
#> [1] 2
```
