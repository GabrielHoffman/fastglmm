# Fitter Function for Linear Mixed Model

Prepare data for model fitting with a call to Rcpp code

## Usage

``` r
fastlmm.fit(
  y,
  X,
  Z,
  offset = NULL,
  REML = FALSE,
  delta = NULL,
  rank = ncol(Z),
  weights = NULL,
  delta.range = c(-10, 10),
  tol = 1e-06,
  lambda = 0,
  nthreads = 6
)
```

## Arguments

- y:

  response vector

- X:

  design matrix

- Z:

  sparse matrix of indicators for random effect

- offset:

  offset

- REML:

  logical scalar - Should the estimates be chosen to optimize the REML
  criterion vs ML?

- delta:

  if `NULL` estimate delta, if value is given used this fixed values

- rank:

  rank of random effect. The maximum rank is the number of columns in
  `Z`. A low rank approximation can be useful if the eigen-values
  decrease quickly.

- weights:

  an optional vector of prior weights with a value for each sample. When
  the response has multiple columns, a vector of weight can be reused
  for each respose, or a matrix the same dimension as the responses
  matrix can weight each response separately.

- delta.range:

  min and max values (in log space), of the search space for delta to
  fit the random effect

- tol:

  convergence criterion for the 1D search of the delta space

- lambda:

  ridge shrinkage parameter

- nthreads:

  number of threads

## Value

`U` and `s` values are from the SVD of weighted Z

## Details

Fit a linear mixed model with a single variance component.
