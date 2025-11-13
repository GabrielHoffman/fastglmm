# Refit fastlmm model with new delta

Refit fastlmm model with new delta value

## Usage

``` r
refit(fit, delta)
```

## Arguments

- fit:

  model fit of class `fastlmm`

- delta:

  new value for delta

## Details

Useful for evaluating log-likelihood and multiple values of delta

## Examples

``` r
library(lme4)

fit <- fastlmm(Reaction ~ Days + (1 | Subject), sleepstudy)

summary(fit)
#> Linear mixed model fit by ML ['fastlmm']
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept) 251.4051     9.5062   26.45   <2e-16 ***
#> Days         10.4673     0.8017   13.06   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Variance components:
#>   sigSq_g: 1297
#>   sigSq_e: 954.5
#>   hSq:     57.6 %
#> 

summary(refit(fit, delta=1000))
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'summary': no applicable method for 'refit' applied to an object of class "fastlmm"
```
