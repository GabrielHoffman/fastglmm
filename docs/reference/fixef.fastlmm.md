# Extract Fixed Effects

Extract Fixed Effects

## Usage

``` r
# S3 method for class 'fastlmm'
fixef(object, ...)
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

fixef(fit)
#>     (Intercept)         trtdrug        trtdrug+ I(week > 2)TRUE 
#>       3.4127450      -1.2475583      -0.7545406      -1.6075955 
```
