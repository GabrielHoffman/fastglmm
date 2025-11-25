# Test Linear Hypothesis

Test Linear Hypothesis

## Usage

``` r
# S3 method for class 'fastlmm'
linearHypothesis(model, ...)

linearHypothesis(model, ...)
```

## Arguments

- model:

  fitted model of class `fastlmm`

- ...:

  other args passed to
  [`car::linearHypothesis.default()`](https://rdrr.io/pkg/car/man/linearHypothesis.html)

## Value

[`car::linearHypothesis()`](https://rdrr.io/pkg/car/man/linearHypothesis.html)

## See also

[`car::linearHypothesis()`](https://rdrr.io/pkg/car/man/linearHypothesis.html)

## Examples

``` r
library(MASS)

# GLMM via PQL
fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

linearHypothesis(fit, "trtdrug", test="F")
#> 
#> Linear hypothesis test:
#> trtdrug = 0
#> 
#> Model 1: restricted model
#> Model 2: y ~ trt + I(week > 2) + (1 | ID)
#> 
#>   Res.Df Df      F  Pr(>F)  
#> 1 190.18                    
#> 2 189.18  1 3.8166 0.05222 .
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
