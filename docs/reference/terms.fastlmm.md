# Model Terms

Model Terms

## Usage

``` r
# S3 method for class 'fastlmm'
terms(x, ...)
```

## Arguments

- x:

  fitted model of class `fastlmm`

- ...:

  other args, not used

## Examples

``` r
library(MASS)

# GLMM via PQL
fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

terms(fit)
#> y ~ trt + I(week > 2)
#> attr(,"variables")
#> list(y, trt, I(week > 2))
#> attr(,"factors")
#>             trt I(week > 2)
#> y             0           0
#> trt           1           0
#> I(week > 2)   0           1
#> attr(,"term.labels")
#> [1] "trt"         "I(week > 2)"
#> attr(,"order")
#> [1] 1 1
#> attr(,"intercept")
#> [1] 1
#> attr(,"response")
#> [1] 1
#> attr(,".Environment")
#> <environment: 0x14cca41b0>
```
