# Extract the modes of the random effect

Extract the conditional modes of the random effect from the model fit

## Usage

``` r
# S3 method for class 'fastlmm'
ranef(object, ...)
```

## Arguments

- object:

  fitted model of class `fastlmm`

- ...:

  other args, not used

## Value

BLUP, modes of the random effect

## Examples

``` r
library(MASS)

# GLMM via PQL
fit <- fastglmm(y ~ trt + I(week > 2) + (1 | ID),
   family = binomial(), data = bacteria)

ranef(fit)
#> $ID
#>           [,1]
#> X01  0.6768634
#> X02 -0.4371422
#> X03  1.4720466
#> X04  0.8295863
#> X05  0.8295863
#> X06  1.2599034
#> X07 -1.7373346
#> X08 -0.6660902
#> X09  0.8295863
#> X10  0.6802351
#> X11  0.8295863
#> X12 -0.5904345
#> X13  0.8295863
#> X14 -1.4169262
#> X15  1.4720466
#> X16 -0.1085555
#> X17 -0.6660902
#> X18 -1.3434724
#> X19 -1.5889409
#> X20  1.0125179
#> X21  0.8295863
#> Y01  0.8295863
#> Y02 -2.3565175
#> Y03 -0.9802696
#> Y04 -1.9621001
#> Y05 -2.6413720
#> Y06  0.6768634
#> Y07 -0.5904345
#> Y08 -0.1821759
#> Y09  0.8008880
#> Y10  0.2551234
#> Y11  0.8295863
#> Y12  0.6802351
#> Y13 -1.9621001
#> Y14  0.4778158
#> Z01  0.9642264
#> Z02  0.8295863
#> Z03  1.2028040
#> Z05 -1.5889409
#> Z06 -1.1466893
#> Z07 -0.5904345
#> Z09  1.0125179
#> Z10  0.8295863
#> Z11  1.4720466
#> Z14 -0.1085555
#> Z15  1.2039370
#> Z19  1.2028040
#> Z20 -1.0585337
#> Z24 -0.9802696
#> Z26 -0.1085555
#> 
```
