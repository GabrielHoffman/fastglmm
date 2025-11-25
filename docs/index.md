  

## Efficiently Fit Generalized Linear Mixed Model with a Single Random Effect on \>1 Million Samples

Fitting linear mixed models (LMMs) and generalized linear mixed models
(GLMMs) on large-scale datasets can be very computationally expensive.
The **fastglmm** package provides a very fast method for fitting LMMs
and GLMMs with a single random effect. This can be 10-1000x faster than
other methods for large datasets.

The **fastglmm** package is designed to be intuitive for users of other
mixed model software such as
[lme4](https://cran.r-project.org/package=lme4), and uses familar mixed
model syntax and R generics. Since **fastglmm** is designed for large
datasets and GLMMs in mind, it uses maximum likelihood rather than REML
to fit LMMs.

### Example Usage

As an introduction, consider fitting a LMM to the `sleepstudy` data
distributed with `lme4`. The model syntax and parameter estimates are
the same as using `lme4::lmer(..., REML=FALSE)`.

``` r

library(fastglmm)
fit <- fastlmm(Reaction ~ Days + (1 | Subject), sleepstudy)

fit
#> 
#> Call:
#> fastlmm(formula = Reaction ~ Days + (1 | Subject), data = sleepstudy)
#> 
#> Coefficients:
#> (Intercept)         Days  
#>      251.41        10.47  
#> 

summary(fit)
#> Linear mixed model fit by ML ['fastlmm']
#>  Formula: Reaction ~ Days + (1 | Subject)
#> 
#> Coefficients:
#>             Estimate Std. Error t value Pr(>|t|)    
#> (Intercept) 251.4051     9.5062   26.45   <2e-16 ***
#> Days         10.4673     0.8017   13.06   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual df: 162.2 
#> 
#> Variance components:
#>   sigSq_g: 1297
#>   sigSq_e: 954.5
#>   delta:   0.736
```

### Installation

``` r

devtools::install_github("GabrielHoffman/fastglmm")
```

## Implementation details

- Analysis of a given response or features is performed on a single
  thread. Analysis is parallelized across responses/features using
  oneAPI Threading Building Blocks
  ([oneTBB](https://uxlfoundation.github.io/oneTBB/))
- Linear algebra is performed with
  [Armadillo](https://arma.sourceforge.net)
- Use of dense and sparse matricies is supported at the C++ level using
  templates
- Interface between R and C++ is implemented in
  [Rcpp](https://www.rcpp.org) and also depends on
  [RcppArmadillo](https://cran.r-project.org/package=RcppArmadillo) and
  [RcppParallel](https://cran.r-project.org/package=RcppParallel)
- Implemented in C++17

### See header-only C++ library [documentation](http://gabrielhoffman.github.io/fastglmm/doxygen/html/index.md)
