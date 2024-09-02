
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
