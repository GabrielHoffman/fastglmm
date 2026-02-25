# Reconstruct indicator matrix from eigen decomp

Reconstruct indicator matrix from eigen decomp

## Usage

``` r
reconstruct_indicator(dcmp, weights = NULL)
```

## Arguments

- dcmp:

  eigen decomp from
  [`indicator_decomp()`](http://gabrielhoffman.github.io/fastglmm/reference/indicator_decomp.md)

- weights:

  vector of weights with a value for each sample. If ommited, weights
  are set to 1.

## Examples

``` r
ID <- factor(sample(LETTERS[1:4], 100, replace=TRUE))
w <- seq(length(ID))
Z <- preprocess_indicator(ID)
dcmp <- indicator_decomp(ID, w)

Z_recon <- reconstruct_indicator( dcmp, w)

range(Z_recon - Z)
#> [1] -1.110223e-16  2.220446e-16
```
