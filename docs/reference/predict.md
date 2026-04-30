# Model Predictions

Model Predictions

## Usage

``` r
# S3 method for class 'fastlmm'
predict(object, newdata = NULL, type = c("response", "terms"), ...)

# S3 method for class 'fastglmm'
predict(object, newdata = NULL, type = c("link", "response", "terms"), ...)
```

## Arguments

- object:

  fitted model of class `fastlmm`

- newdata:

  optionally, a data.frame in which to look for variables with which to
  predict. If omitted, the fitted predictors are used.

- type:

  the type of prediction required. The default is on the scale of the
  linear predictors; the alternative `"response"` is on the scale of the
  response variable. Thus for a default binomial model the default
  predictions are of log-odds (probabilities on logit scale) and
  `type = "response"` gives the predicted probabilities.
  `type = "terms"` computes the linear predict for each model term

- ...:

  other args, not used

## Examples

``` r
library(MASS)

# GLMM via PQL
fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

predict(fit)[1:3]
#> [1] 4.048960 4.048960 2.441364
```
