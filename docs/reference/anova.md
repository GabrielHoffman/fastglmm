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

## Value

anova table like [`stats::anova()`](https://rdrr.io/r/stats/anova.html)

## Examples

``` r
library(MASS)

# GLMM via PQL
fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

anova(fit)
#> Analysis of Variance Table
#>             df1    df2      F   Pr(>F)    
#> (Intercept)   1 223.72 25.791 8.01e-07 ***
#> trt           2 108.33  1.111 0.332967    
#> I(week > 2)   1 490.21 12.655 0.000411 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
