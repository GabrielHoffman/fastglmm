# ANOVA Tables

ANOVA Tables

## Usage

``` r
# S3 method for class 'fastlmm'
anova(object, ddf = c("satterthwaite", "asymptotic"), ...)
```

## Arguments

- object:

  fitted model of class `fastlmm`

- ddf:

  `"satterthwaite"`: use Satterthwaite approximation to denominator
  degrees of freedom for the F distribution, or `"asymptotic"` to use
  chisq distribution as null for the test statistic

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
#>             df1    df2      F    Pr(>F)    
#> (Intercept)   1 267.08 27.123 3.817e-07 ***
#> trt           2 135.49  2.067 0.1305425    
#> I(week > 2)   1 485.74 12.609 0.0004213 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
