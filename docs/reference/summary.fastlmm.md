# Object Summaries and Hypothesis Testing

Object summaries and hypothesis testing of fixed effects

## Usage

``` r
# S3 method for class 'fastlmm'
summary(object, ddf = c("satterthwaite", "asymptotic"), ...)
```

## Arguments

- object:

  fitted model of class `fastlmm`

- ddf:

  `"satterthwaite"`: use Satterthwaite approximation to denominator
  degrees of freedom for the Student-t distribution, or `"asymptotic"`
  to use normal distribution as null for the test statistic

- ...:

  other args, not used

## Value

summary of model fit

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
```
