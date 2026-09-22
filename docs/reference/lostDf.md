# Degrees of freedom lost

Degrees of freedom lost

## Usage

``` r
lostDf(fit, cutoff = 1e-04)
```

## Arguments

- fit:

  model fit

- cutoff:

  cutoff for fitted values

## Value

degrees of freedom lost due to 0 counts in response

## References

Lun, A. T., & Smyth, G. K. (2017). No counts, no variance: allowing for
loss of degrees of freedom when assessing biological variability from
RNA-seq data. Statistical Applications in Genetics and Molecular
Biology, 16(2), 83-93.

## Examples

``` r
library(MASS)

# GLMM via PQL
fit = fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

fastglmm:::lostDf(fit)
#> [1] 0
```
