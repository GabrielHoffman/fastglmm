
<br>

## Efficient Fit of Linear Mixed Model with a Single Random Effect on >1 Million Samples

Fitting linear mixed model on large-scale datasets can be very computationally expensive.  Here, we provide a very fast method for fitting linear mixed models with a single random effect.  This can be 10-100x faster than other methods for large datasets.



```r
devtools::install_github("GabrielHoffman/fastglmm")
```

### See header-only C++ library [documentation](doxygen/html/index.html)

## Implmentation details
- Analysis of a given response or features is performed on a single thread.  Analysis is parallelized across responses/features using oneAPI Threading Building Blocks ([oneTBB](https://uxlfoundation.github.io/oneTBB/))
- Linear algebra is performed with [Armadillo](https://arma.sourceforge.net)
- Use of dense and sparse matricies is supported at the C++ level uses C++ templates
- Interface between R and C++ is implemented in [Rcpp](https://www.rcpp.org) and also depends on [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/) and [RcppParallel](https://cran.r-project.org/web/packages/RcppParallel/)
- Implemented in C++17