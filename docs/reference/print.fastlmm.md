# Print model

Print model

## Usage

``` r
# S3 method for class 'fastlmm'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  fitted model of class `fastlmm`

- digits:

  minimal number of \_significant\_ digits,

- ...:

  other args, not used

## Value

print model

## Examples

``` r
data(PsychAD)

# regression formula
form <- PTPRG ~ (1|SubID) + offset(log(libSize))

# NB GLMM on PTPRG expression via PQL
fit <- fastglmm.nb(form, PsychAD)

fit
#> 
#> Call:
#> fastglmm(formula = formula, data = data, family = negative.binomial(NA), 
#>     weights = weights, maxit = maxit, tol = tol, tol.eta = tol.eta, 
#>     doCoxReid = doCoxReid, lambda = lambda, nthreads = nthreads)
#> 
#> Coefficients:
#> (Intercept)  
#>      -8.639  
#> 
```
