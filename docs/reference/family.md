# Extract Model Family

Extract Model Family

## Usage

``` r
# S3 method for class 'fastlmm'
family(object, ...)

# S3 method for class 'fastglmm'
family(object, ...)

# S3 method for class 'glmmPQL'
family(object, ...)
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

family(fit)
#> 
#> Family: binomial 
#> Link function: logit 
#> 
```
