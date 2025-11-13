# Changelog

## fastglmm 0.3.3

- Oct 17, 2025
- when model fit fails return NaN values
- Update docs

## fastglmm 0.3.2

- Oct 7, 2025
- add checks and compatibility with BatchRegression

## fastglmm 0.3.1

- Sept 16, 2025
- add and check generics
- additional testing

## fastglmm 0.3.0

- May 6, 2025
- rename
- add
  [`fastglmm()`](http://gabrielhoffman.github.io/fastglmm/reference/fastglmm.md)
  in C++
- pull code from `BatchRegression`

## fastglmm 0.2.1

- April 24, 2025
- refactor for compatibility across ecosystem

## fastglmm 0.2.0

- March 6, 2025
- Major refactor

## fastglmm 0.1.6

- Nov 19, 2024
- fix API for linear and mixed model regression

## fastglmm 0.1.5

- Nov 7, 2024
- `linearRegression.h` supports weighted regression with preprojection

## fastglmm 0.1.4

- Sept 10, 2024
- move `fastlmmLib` code to `inst/include` for accessible header-only
  library

## fastglmm 0.1.3

- Aug 28, 2024
- changes to allow PQL with
  [`fastglmm()`](http://gabrielhoffman.github.io/fastglmm/reference/fastglmm.md)

## fastglmm 0.1.2

- Aug 7, 2024
- multivariate model is run in parallel

## fastglmm 0.1.1

- Aug 6, 2024
- fix inconsistent merge

## fastglmm 0.1.0

- Aug 6, 2024
- fix bottleneck in
  - [`model.frame()`](https://rdrr.io/r/stats/model.frame.html) for
    matrix response
  - log-likelihood in Rcpp since `weights` are constant across
    iterations
  - memory usage
- pass R CMD check
- RcppArmadillo code works except doesn’t consider varying weights
  across responses
