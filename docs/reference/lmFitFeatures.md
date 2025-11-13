# Fit series of linear regression models with the same response

Fit regression model `y ~ design + X_features[,j]` for each feature j

## Usage

``` r
lmFitFeatures(
  y,
  design,
  data,
  weights,
  detail = 0,
  preprojection = TRUE,
  nthreads = 1,
  ...
)

# S4 method for class 'ANY,ANY,matrix'
lmFitFeatures(
  y,
  design,
  data,
  weights,
  detail = 0,
  preprojection = TRUE,
  nthreads = 1,
  ...
)
```

## Arguments

- y:

  response vector

- design:

  design matrix shared across all models

- data:

  feature matrix with model j using feature j

- weights:

  sample-level weights

- detail:

  level of model detail returned, with LOW = 0, MEDIUM = 1, HIGH = 2.
  LOW (beta, se, sigSq, rdf), MEDIUM (vcov), HIGH (residuals), MOST
  (hatvalues)

- preprojection:

  default TRUE. Use preproject of design matrix to accelerate
  calculations

- nthreads:

  number of threads. Each model is fit in serial, analysis is
  parallelized across features

- ...:

  other args

## Value

List of parameter estimates with entries `coef`, `se`, `sigSq`, `rdf`
and other depending on `detail`

## Examples

``` r
n <- 100 # number of samples
p <- 10 # number of features
nc <- 3 # number shared covariates
set.seed(1)
y <- rnorm(n)
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- seq(p)
design <- matrix(rnorm(n * nc), n, nc)
w <- seq(n)
w <- w / mean(w)

# fit regressions with model j including X[,j]
fit <- lmFitFeatures(y, design, X, w)

# examine results
lapply(fit, head, 2)
#> $coef
#>          [,1]
#> 1 -0.01673020
#> 2  0.02287094
#> 
#> $se
#>         [,1]
#> 1 0.09250378
#> 2 0.09378003
#> 
#> $sigSq
#>         1         2 
#> 0.8397724 0.8395384 
#> 
#> $rdf
#>  1  2 
#> 96 96 
#> 
```
