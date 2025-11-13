# Hypothesis test for heritability

Hypothesis test for heritability using either asymptotic normal
approximation or empirical permutations.

## Usage

``` r
heritability(fit, method = c("information", "permutation"), nperms = 100)
```

## Arguments

- fit:

  model fit of class `fastlmm`

- method:

  `"information"` or `"permutation"`,

- nperms:

  number of permutations

## Details

For `method == "information"`, the profile log-likelihood is evaluted
with respect to hsq. The standard error is obtained from the information
matrix based on the Hessian evaluated at the MLE of hsq. The p-value is
then computed from this estimated standard error using a normal
approximation. This approach can perform well for large sample sizes, by
performs poorly for moderate sample sizes.

For `method == "permutation"`, ....

## References

Abney, M. (2015). Permutation testing in the presence of polygenic
variation. Genetic epidemiology, 39(4), 249-258.
<https://doi.org/10.1002/gepi.21893>
