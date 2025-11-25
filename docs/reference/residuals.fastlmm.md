# Extract Model Residuals

Extract Model Residuals

## Usage

``` r
# S3 method for class 'fastlmm'
residuals(object, type = c("working", "response", "deviance", "pearson"), ...)
```

## Arguments

- object:

  fitted model of class `fastlmm`

- type:

  the type of residuals which should be returned.

- ...:

  other args, not used

## See also

[`residuals.glm()`](https://rdrr.io/r/stats/glm.summaries.html)

## Examples

``` r
library(MASS)

# GLMM via PQL
fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

residuals(fit)[1:3]
#> [1] 0.1859579 0.1859579 0.4085591
```
