# Denominator Degrees of Freedom

Denominator degrees of freedom using Satterthwaite method

## Usage

``` r
ddf(fit, L = diag(1, length(coef(fit))))
```

## Arguments

- fit:

  model fit from
  [`fastlmm()`](http://gabrielhoffman.github.io/fastglmm/reference/fastlmm.md)
  or
  [`fastglmm()`](http://gabrielhoffman.github.io/fastglmm/reference/fastglmm.md)

- L:

  matrix of coefficient contrasts, one per \_row\_

## Value

array, denominator degrees of freedom for each contrast (i.e. \_row\_)

## Examples

``` r
library(MASS)

# GLMM via PQL
fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

# denominator degrees of freedom 
# used for hypothesis testing below
ddf(fit)
#> [1] 267.0785 123.8255 131.7683 485.7446

# summarize fit
summary(fit)
#> Generalized linear mixed model fit by PQL ['fastglmm']
#>  Family: binomial  ( logit )
#>  Formula: y ~ trt + I(week > 2) + (1 | ID)
#> 
#> Coefficients:
#>                 Estimate Std. Error   ddf t value Pr(>|t|)    
#> (Intercept)       3.4127     0.6553 267.1   5.208 3.82e-07 ***
#> trtdrug          -1.2476     0.8141 123.8  -1.533 0.127950    
#> trtdrug+         -0.7545     0.8157 131.8  -0.925 0.356663    
#> I(week > 2)TRUE  -1.6076     0.4527 485.7  -3.551 0.000421 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual df: 189.2 
#> 
#> Variance components:
#>   sigSq_g: 1.993
#>   sigSq_e: 5.117
#>   delta:   2.568
```
