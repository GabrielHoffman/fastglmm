# Test Linear Hypothesis

Test Linear Hypothesis

## Usage

``` r
linearHypothesis(model, ...)

# S3 method for class 'fastlmm'
linearHypothesis(
  model,
  hypothesis.matrix,
  rhs = NULL,
  ...,
  ddf = c("satterthwaite", "asymptotic")
)
```

## Arguments

- model:

  fitted model of class `fastlmm`

- ...:

  other args passed to
  [`car::linearHypothesis.default()`](https://rdrr.io/pkg/car/man/linearHypothesis.html)

- hypothesis.matrix:

  matrix (or vector) giving linear combinations of coefficients by rows,
  or a character vector giving the hypothesis in symbolic form

- rhs:

  right-hand-side vector for hypothesis, with as many entries as rows in
  the hypothesis matrix; can be omitted, in which case it defaults to a
  vector of zeroes. For a multivariate linear model, ‘rhs’ is a matrix,
  defaulting to 0

- ddf:

  `"satterthwaite"`: use Satterthwaite approximation to denominator
  degrees of freedom for the Student-t or F distribution, or
  `"asymptotic"` to use normal distribution or chisq as null for the
  test statistic

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
#>   Res.Df Df      F Pr(>F)
#> 1 124.83                 
#> 2 123.83  1 2.3486 0.1279
```
