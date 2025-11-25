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
#> (Intercept)  1 44.078  3.156e-11 ***
#> trt          2  3.359     0.1865    
#> I(week > 2)  1 20.491  5.993e-06 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
