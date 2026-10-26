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

## Value

refit model

## Details

Useful for evaluating log-likelihood and multiple values of delta

## Examples

``` r
library(lme4)

fit <- fastlmm(Reaction ~ Days + (1 | Subject), sleepstudy)

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

summary(refitModel(fit, delta=1000))
#> Linear mixed model fit by REML ['fastlmm']
#>  Formula: Reaction ~ Days + (1 | Subject)
#> 
#> Coefficients:
#>             Estimate Std. Error    df t value Pr(>|t|)    
#> (Intercept)  251.405      6.599 179.8   38.10  < 2e-16 ***
#> Days          10.467      1.234 178.7    8.48 8.29e-15 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual df: 177.8 
#> 
#> Variance components:
#>   sigSq_g: 2.263
#>   sigSq_e: 2263
#>   delta:   1000
```
