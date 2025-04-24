
# fastlmm 0.2.1
 - April 24, 2025
 - refactor for compatibility across ecosystem

# fastlmm 0.2.0
 - March 6, 2025
 - Major refactor

# fastlmm 0.1.6
 - Nov 19, 2024
 - fix API for linear and mixed model regression

# fastlmm 0.1.5
 - Nov 7, 2024
 - `linearRegression.h` supports weighted regression with preprojection


# fastlmm 0.1.4
 - Sept 10, 2024
 - move `fastlmmLib` code to `inst/include` for accessible header-only library

# fastlmm 0.1.3
 - Aug 28, 2024
 - changes to allow PQL with `fastglmm()`

# fastlmm 0.1.2
 - Aug 7, 2024
 - multivariate model is run in parallel

# fastlmm 0.1.1
 - Aug 6, 2024
 - fix inconsistent merge

# fastlmm 0.1.0
 - Aug 6, 2024
 - fix bottleneck in 
 	- `model.frame()` for matrix response
 	- log-likelihood in Rcpp since `weights` are constant across iterations
 	- memory usage
 - pass R CMD check
 - RcppArmadillo code works except doesn't consider varying weights across responses
