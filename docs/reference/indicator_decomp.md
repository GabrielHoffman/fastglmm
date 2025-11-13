# Spectral decomposition of factor indicator matrix

Given a factor, construct the spectral decompostion of the corresponding
indicator matrix. Uses linear time algorithm and sparse matrix algebra.
Resulting vector space is also sparse.

## Usage

``` r
indicator_decomp(x, weights = NULL, rank = NULL, sort = FALSE)
```

## Arguments

- x:

  object of type `factor` or `sparseMatrix`

- weights:

  vector of weights with a value for each sample. If ommited, weights
  are set to 1.

- rank:

  rank of random effect. The maximum rank is the number of columns in
  `Z`. A low rank approximation can be useful if the eigen-values
  decrease quickly.

- sort:

  sort eigen vectors and values

## Value

list storing spectral decomposition like from
[`eigen()`](https://rdrr.io/r/base/eigen.html)

## Details

This approach is dramatically faster than the naive algorithm that is
quadratic time in the number of levels.

## See also

[`preprocess_indicator()`](http://gabrielhoffman.github.io/fastglmm/reference/preprocess_indicator.md)

## Examples

``` r
set.seed(2)

# get 10 samples from 4 individuals
indiv <- factor(sample(LETTERS[seq(4)], 10, replace = TRUE))

# fast spectral decomposition
indicator_decomp(indiv, sort = TRUE)
#> $vectors
#> 10 x 4 sparse Matrix of class "dgCMatrix"
#>         A         D         B C
#>  [1,] 0.5 .         .         .
#>  [2,] .   .         .         1
#>  [3,] .   .         0.7071068 .
#>  [4,] .   .         0.7071068 .
#>  [5,] .   0.5773503 .         .
#>  [6,] .   0.5773503 .         .
#>  [7,] 0.5 .         .         .
#>  [8,] 0.5 .         .         .
#>  [9,] 0.5 .         .         .
#> [10,] .   0.5773503 .         .
#> 
#> $values
#> A D B C 
#> 4 3 2 1 
#> 

# Create spectral decomposition
# using a slower approach ignoring the
# special structure of the problem.
# This is quadratic time in the number of levels
library(Matrix)
Z <- t(fac2sparse(indiv))
ZtZ <- crossprod(Z)
dcmp <- eigen(ZtZ, symmetric = TRUE)
dcmp$vectors <- Matrix(zapsmall(dcmp$vectors), sparse = TRUE)
A <- solve(dcmp$vectors)
D <- Diagonal(length(dcmp$values), 1 / sqrt(dcmp$values))
dcmp$vectors <- tcrossprod(Z, A) %*% D

dcmp
#> eigen() decomposition
#> $values
#> [1] 4 3 2 1
#> 
#> $vectors
#> 10 x 4 sparse Matrix of class "dgCMatrix"
#>                                
#>  [1,] 0.5 .         .         .
#>  [2,] .   .         .         1
#>  [3,] .   .         0.7071068 .
#>  [4,] .   .         0.7071068 .
#>  [5,] .   0.5773503 .         .
#>  [6,] .   0.5773503 .         .
#>  [7,] 0.5 .         .         .
#>  [8,] 0.5 .         .         .
#>  [9,] 0.5 .         .         .
#> [10,] .   0.5773503 .         .
#> 
```
