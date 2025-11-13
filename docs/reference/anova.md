# ANOVA Tables

ANOVA Tables

## Usage

``` r
# S3 method for class 'fastlmm'
anova(object, ...)

# S3 method for class 'fastglmm'
anova(object, ...)
```

## Arguments

- object:

  fitted model of class `fastlmm`

- ...:

  other args, not used

## Examples

``` r
library(MASS)

# GLMM via PQL
fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

anova(fit)
#> Analysis of Variance Table
#>             df  Chisq Pr(>Chisq)    
#> (Intercept)  1 44.084  3.145e-11 ***
#> trt          2  3.360     0.1864    
#> I(week > 2)  1 20.486  6.007e-06 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
