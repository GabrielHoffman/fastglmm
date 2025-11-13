# Diagnostic Plots for Model Fits

Diagnostic plots for `fastlmm` Fits

## Usage

``` r
# S3 method for class 'fastlmm'
plot(x, form = resid(x, type = "pearson") ~ fitted(x), ...)
```

## Arguments

- x:

  fitted model of class `fastlmm`

- form:

  formula to plot

- ...:

  other args, not used

## Examples

``` r
library(MASS)

# GLMM via PQL
fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

plot(fit)
```
