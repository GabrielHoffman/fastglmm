# Refit fastlmm model with new delta

Refit fastlmm model with new delta value

## Usage

``` r
refitModel(fit, delta = NULL, interceptOnly = FALSE, fixedNBtheta = FALSE)
```

## Arguments

- fit:

  model fit of class `fastlmm` or `fastglmm`

- delta:

  new value for delta

- interceptOnly:

  if `TRUE`, only retain the intercept from the fixed effect design
  matrix

- fixedNBtheta:

  if `FALSE`, low theta in NB model to be re-estimated too

## Details

Useful for evaluating log-likelihood and multiple values of delta

## Examples

``` r
library(lme4)

fit <- fastlmm(Reaction ~ Days + (1 | Subject), sleepstudy)

summary(fit)
#> Linear mixed model fit by ML ['fastlmm']
#>  Formula: Reaction ~ Days + (1 | Subject)
#> 
#> Coefficients:
#>             Estimate Std. Error    df t value Pr(>|t|)    
#> (Intercept) 251.4051     9.5062  24.5   26.45   <2e-16 ***
#> Days         10.4673     0.8017 162.0   13.06   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual df: 162.2 
#> 
#> Variance components:
#>   sigSq_g: 1297
#>   sigSq_e: 954.5
#>   delta:   0.736

summary(refitModel(fit, delta=1000))
#> Linear mixed model fit by ML ['fastlmm']
#>  Formula: Reaction ~ Days + (1 | Subject)
#> 
#> Coefficients:
#>             Estimate Std. Error    df t value Pr(>|t|)    
#> (Intercept)  251.405      6.563 127.9  38.308  < 2e-16 ***
#> Days          10.467      1.228 162.0   8.527 1.01e-14 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual df: 177.8 
#> 
#> Variance components:
#>   sigSq_g: 2.238
#>   sigSq_e: 2238
#>   delta:   1000
```
