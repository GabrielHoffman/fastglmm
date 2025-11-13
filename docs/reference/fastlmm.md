# Fast Linear Mixed Model with 1 Random eEffect

Fit a linear mixed-effects model with 1 random effect with REML or
maximum likelihood using the very fast algorithm and implementation

Stores results of fastlmm model fit

## Usage

``` r
fastlmm(
  formula,
  data,
  REML = FALSE,
  delta = NULL,
  weights = NULL,
  delta.range = c(-10, 10),
  tol = 1e-06,
  nthreads = 6
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

- REML:

  logical scalar - Should the estimates be chosen to optimize the REML
  criterion vs ML?

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

## Value

none

## Details

Hoffman (2013), Lippert, et al. (2011)

## References

Hoffman, G. E. (2013). Correcting for population structure and kinship
using the linear mixed model: theory and extensions. PloS one, 8(10),
e75707. <https://doi.org/10.1371/journal.pone.0075707>

Lippert, C., Listgarten, J., Liu, Y., Kadie, C. M., Davidson, R. I., &
Heckerman, D. (2011). FaST linear mixed models for genome-wide
association studies. Nature methods, 8(10), 833-835.
<https://www.nature.com/articles/nmeth.1681>

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
#> Linear mixed model fit by ML ['fastlmm']
#>  Formula: Reaction ~ Days + (1 | Subject)
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept) 251.4051     9.5062   26.45   <2e-16 ***
#> Days         10.4673     0.8017   13.06   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual df: 162.2 
#> 
#> Variance components:
#>   sigSq_g: 1297
#>   sigSq_e: 954.5
#>   delta:   0.736

# results are identical for lmer(...,REML=FALSE)
fit2 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy, REML = FALSE)
coef(summary(fit2))
#>              Estimate Std. Error  t value
#> (Intercept) 251.40510  9.5061852 26.44648
#> Days         10.46729  0.8017354 13.05579
```
