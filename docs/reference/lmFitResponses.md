# Fit series of linear regression models to multiple responses with shared design matrix

Fit regression model `Y[j,] ~ design` for each feature j

## Usage

``` r
lmFitResponses(Y, design, Weights, detail = 0, nthreads = 1, ...)

# S4 method for class 'matrix'
lmFitResponses(Y, design, Weights, detail = 0, nthreads = 1, ...)
```

## Arguments

- Y:

  matrix of responses as \_\_rows\_\_

- design:

  design matrix

- Weights:

  matrix sample-level weights the same dimension as Y

- detail:

  level of model detail returned, with LOW = 0, MEDIUM = 1, HIGH = 2.
  LOW (`beta`, `se`, `sigSq`, `rdf`), MEDIUM (`vcov`), HIGH
  (`residuals`), MOST (`hatvalues`)

- nthreads:

  number of threads. Each model is fit in serial, analysis is
  parallelized across responses.

- ...:

  other args

## Value

List of parameter estimates with entries `coef`, `se`, `sigSq`, `rdf`
and other depending on `detail`

## Details

Since the weights vary for each response, each model is computed
separately without recycling precomputed values

## Examples

``` r
n <- 100
m <- 5
nc <- 2
set.seed(1)
Y <- matrix(rnorm(n * m), m, n)
X <- matrix(rnorm(n * nc), n, nc)
rownames(Y) <- seq(m)
W <- matrix(runif(n * m), m, n)

# fit regressions with model j using Y[,j] as a response
fit <- lmFitResponses(Y, X, W)

# examine results
lapply(fit, head, 2)
#> $coef
#>         [,1]        [,2]
#> 1 -0.2167269 -0.26282250
#> 2  0.1217050  0.08171999
#> 
#> $se
#>         [,1]       [,2]
#> 1 0.09798615 0.08866582
#> 2 0.09649571 0.08894915
#> 
#> $sigSq
#>         1         2 
#> 0.4200639 0.4111836 
#> 
#> $rdf
#>  1  2 
#> 98 98 
#> 
```
