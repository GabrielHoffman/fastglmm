# Fit linear mixed model using SVD of covariance

Fit linear mixed model using SVD of covariance to scale to large sample
sizes.

## Usage

``` r
fastlmm3(
  Y,
  X,
  U,
  s,
  delta = NULL,
  sig_a_fixed = FALSE,
  rank = ncol(U),
  W_til = NULL
)
```

## Arguments

- Y:

  response vector

- X:

  matrix of covariates

- U:

  principal components of covariance matrix

- s:

  eigen values from of covariance matrix

- delta:

  ratio of variance components estimated using

- rank:

  number of of principal components used

- W_til:

  components to remove from U. (Not recommended)

## Value

summary statistics for model fit, and hypothesis testing using
X_test_lst, if available

## Details

Fit a linear mixed model with a single variance component.
