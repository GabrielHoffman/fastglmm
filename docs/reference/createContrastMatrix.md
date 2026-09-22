# Create Contrast Matrix

Create contrast matrix as array or text

## Usage

``` r
createContrastMatrix(model, hypothesis.matrix, rhs = NULL)
```

## Arguments

- model:

  model fit

- hypothesis.matrix:

  array or text indicating contrast

- rhs:

  right hand side of equation

## Value

contrast matrix

## See also

[`car::linearHypothesis.default()`](https://rdrr.io/pkg/car/man/linearHypothesis.html)

## Examples

``` r
library(lme4)
#> Loading required package: Matrix
#> 
#> Attaching package: ‘lme4’
#> The following objects are masked from ‘package:fastglmm’:
#> 
#>     getLambda, getTheta
#> The following object is masked from ‘package:nlme’:
#> 
#>     lmList

fit <- fastlmm(Reaction ~ Days + (1 | Subject), sleepstudy)

createContrastMatrix(fit, "Days = 0")
#>          (Intercept) Days
#> Days = 0           0    1
```
