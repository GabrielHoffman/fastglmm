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
#> [1] 223.7220 105.3093 111.5378 490.2091

# summarize fit
summary(fit)
#> Generalized linear mixed model fit by PQL ['fastglmm']
#>  Family: binomial  ( logit )
#>  Formula: y ~ trt + I(week > 2) + (1 | ID)
#> 
#> Coefficients:
#>                 Estimate Std. Error    df t value Pr(>|t|)    
#> (Intercept)       3.4756     0.6844 223.7  5.0785 8.01e-07 ***
#> trtdrug          -1.2652     0.8610 105.3 -1.4693 0.144726    
#> trtdrug+         -0.7730     0.8612 111.5 -0.8976 0.371328    
#> I(week > 2)TRUE  -1.6363     0.4600 490.2 -3.5574 0.000411 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual df: 188.2 
#> 
#> Variance components:
#>   sigSq_g: 2.29
#>   sigSq_e: 5.22
#>   delta:   2.28
```
