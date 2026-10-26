# Fast Linear Mixed Model with 1 Random Effect

Fit a linear mixed-effects model with 1 random effect with maximum
likelihood or REML using the very fast algorithm and implementation

## Usage

``` r
fastlmm(
  formula,
  data,
  delta = NULL,
  weights = NULL,
  delta.range = c(-10, 10),
  tol = 1e-06,
  nthreads = 6,
  lambda = 0,
  REML = TRUE
)
```

## Arguments

- formula:

  a two-sided linear formula object describing both the fixed-effects
  and random-effects part of the model, with the response on the left of
  a `~` operator and the terms, separated by `+` operators, on the
  right. Random-effects terms are distinguished by vertical bars (`|`)
  separating expressions for design matrices from grouping factors.

- data:

  an optional data frame containing the variables named in

- delta:

  if `NULL` estimate delta, if value is given use this fixed value

- weights:

  an optional vector of prior weights with a value for each sample. When
  the response has multiple columns, a vector of weight can be reused
  for each respose, or a matrix the same dimension as the responses
  matrix can weight each response separately.

- delta.range:

  min and max values (in log space), of the search space for delta to
  fit the random effect

- tol:

  convergence criterion for the 1D search of the delta space

- nthreads:

  number of threads

- lambda:

  ridge shrinkage parameter

- REML:

  use REML or ML estimates (default to `TRUE`)

## Value

model fit of class fastlmm

## Details

See Hoffman (2013) and Lippert, et al. (2011) for algorithm for low rank
linear mixed model with one random effect

## References

Hoffman, G. E. (2013). Correcting for population structure and kinship
using the linear mixed model: theory and extensions. PloS one, 8(10),
e75707.
[doi:10.1371/journal.pone.0075707](https://doi.org/10.1371/journal.pone.0075707)

Lippert, C., Listgarten, J., Liu, Y., Kadie, C. M., Davidson, R. I., &
Heckerman, D. (2011). FaST linear mixed models for genome-wide
association studies. Nature methods, 8(10), 833-835.
[doi:10.1038/nmeth.1681](https://doi.org/10.1038/nmeth.1681)

## See also

[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html)

## Examples

``` r
library(lme4)

fit <- fastlmm(Reaction ~ Days + (1 | Subject), sleepstudy)

fit
#> 
#> Call:
#> fastlmm(formula = Reaction ~ Days + (1 | Subject), data = sleepstudy)
#> 
#> Coefficients:
#> (Intercept)         Days  
#>      251.41        10.47  
#> 

summary(fit)
#> Linear mixed model fit by REML ['fastlmm']
#>  Formula: Reaction ~ Days + (1 | Subject)
#> 
#> Coefficients:
#>             Estimate Std. Error     df t value Pr(>|t|)    
#> (Intercept) 251.4051     9.7467  22.62   25.79   <2e-16 ***
#> Days         10.4673     0.8042 161.00   13.02   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual df: 162.1 
#> 
#> Variance components:
#>   sigSq_g: 1378
#>   sigSq_e: 960.5
#>   delta:   0.6969

# results are identical for lmer()
fit2 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
coef(summary(fit2))
#>              Estimate Std. Error  t value
#> (Intercept) 251.40510  9.7467163 25.79383
#> Days         10.46729  0.8042214 13.01543
```
