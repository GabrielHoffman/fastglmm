# Fit linear mixed model using SVD of covariance

Fit linear mixed model using SVD of covariance to scale to large sample
sizes.

## Usage

``` r
fastlmm_R(
  Y,
  X,
  U,
  s,
  weights = rep(1, nrow(X)),
  Xu = NULL,
  Yu = NULL,
  delta = NULL,
  sig_a_fixed = FALSE,
  rank = ncol(U)
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

- weights:

  vector weights with value for each sample

- Xu:

  pre-transformed X value

- Yu:

  pre-transformed Y value

- delta:

  ratio of variance components estimated using

- sig_a_fixed:

  if `FALSE`, estimate `sigSq_a` from data

- rank:

  number of of principal components used

## Value

summary statistics for model fit, and hypothesis testing

## Details

Fit a linear mixed model with a single variance component.

## Examples

``` r
library(lme4)
fit <- fastlmm(Reaction ~ Days + (1 | Subject), sleepstudy)

y <- sleepstudy$Reaction
X <- model.matrix(~Days, sleepstudy)
dcmp <- indicator_decomp(sleepstudy$Subject)

fit2 <- fastlmm_R(y, X, dcmp$vectors, dcmp$values)

fit2
#> 
#> Call:
#> NULL
#> 
#> Coefficients:
#> (Intercept)         Days  
#>      251.41        10.47  
#> 
```
