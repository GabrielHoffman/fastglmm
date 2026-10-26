# Calculate Variance-Covariance Matrix for a Fitted Model Object

Calculate Variance-Covariance Matrix for a Fitted Model Object

## Usage

``` r
# S3 method for class 'fastlmm'
vcov(object, ...)
```

## Arguments

- object:

  fitted model of class `fastlmm`

- ...:

  other args, not used

## Value

variance-covariance matrix of fixed effects

## Examples

``` r
library(MASS)

# GLMM via PQL
fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

vcov(fit)
#>                 (Intercept)     trtdrug      trtdrug+ I(week > 2)TRUE
#> (Intercept)       0.4683890 -0.35446875 -0.3398799959   -0.1644577324
#> trtdrug          -0.3544687  0.74140314  0.3404880521    0.0178916075
#> trtdrug+         -0.3398800  0.34048805  0.7417217682   -0.0008778027
#> I(week > 2)TRUE  -0.1644577  0.01789161 -0.0008778027    0.2115859390
```
