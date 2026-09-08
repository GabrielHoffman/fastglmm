# Changelog

## fastglmm 0.4.14

- Sept 6, 2026
- update docs

## fastglmm 0.4.13

- July 28, 2026
- fix accouting for count model offset

## fastglmm 0.4.12

- July 24, 2026
- `_log_moments_nb_mu()` uses approximation for large mu values

## fastglmm 0.4.11

- July 15, 2026
- fix bug in
  [`getTheta()`](http://gabrielhoffman.github.io/fastglmm/reference/getTheta.md)
- use weighted variances in `vpOther()`
- bug fixes

## fastglmm 0.4.10

- July 9, 2026
- [`varpart()`](http://gabrielhoffman.github.io/fastglmm/reference/varpart.md)
  works with [`lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html),
  [`glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html), `glmmTMB()`,
  [`glm()`](https://rdrr.io/r/stats/glm.html),
  [`lm()`](https://rdrr.io/r/stats/lm.html)
  - also
    [`varianceTerms()`](http://gabrielhoffman.github.io/fastglmm/reference/varianceTerms.md),
    `predictTerms()`
- pseudocount defaults to 1
- improve numerical stability of
  [`log_moments_nb_XB()`](http://gabrielhoffman.github.io/fastglmm/reference/log_moments_nb_XB.md)

## fastglmm 0.4.9

- June 17, 2026
- new arguments to
  [`varpart()`](http://gabrielhoffman.github.io/fastglmm/reference/varpart.md)
  - add `log_moments_nb()`

## fastglmm 0.4.8

- June 3, 2026
- add “exact” method for
  [`varpart()`](http://gabrielhoffman.github.io/fastglmm/reference/varpart.md)

## fastglmm 0.4.7

- May 19, 2026
- update count filtering

## fastglmm 0.4.6

- April 24, 2026
- for count regression model, throws error for non-integer responses

## fastglmm 0.4.5

- April 20, 2026
- [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) accepts
  `newdata` parameter

## fastglmm 0.4.4

- April 17, 2026
- fix computing mu (fitted.values) for glmmFitResponses

## fastglmm 0.4.3

- March 31, 2026
- add , ,

## fastglmm 0.4.2

- March 23, 2026
- [`ddf()`](http://gabrielhoffman.github.io/fastglmm/reference/ddf.md)
  now has a min value of 2

## fastglmm 0.4.1

- Feb 26, 2026
- fix hypothesis testing

## fastglmm 0.4.0

- Feb 17, 2026
- Hypothesis test uses Satterthwaite denominator degrees of freedom
  - updated [`summary()`](https://rdrr.io/r/base/summary.html),
    [`anova()`](https://rdrr.io/r/stats/anova.html) and
    [`linearHypothesis()`](http://gabrielhoffman.github.io/fastglmm/reference/linearHypothesis.md)
  - implemented at C++ level, return `A.sat`, `B.sat`, `hessian.vc`

## fastglmm 0.3.9

- Feb 11, 2026
- improved `get_rdf()` for `fastlmm` in C++
- clean up code, consts
- add mutex to response models
- fix bug in
  [`dispersion()`](http://gabrielhoffman.github.io/fastglmm/reference/dispersion.md),
  uses Pearson residuals
- NB-QL
- `negative.binomial(NA)` fixes dispersion to 1
- `negative.binomial(theta)` estimates dispersion from pearson residuals
  and scales vcov and se
- in C++ `NB()` now uses QL dispersion of theta is given

## fastglmm 0.3.8

- Jan 16, 2026
- more flexable
  [`isCountModel()`](http://gabrielhoffman.github.io/fastglmm/reference/isCountModel.md)
  and
  [`isNB()`](http://gabrielhoffman.github.io/fastglmm/reference/isNB.md)
- fix [`simulate()`](https://rdrr.io/r/stats/simulate.html) for NB
  models

## fastglmm 0.3.7

- performance improvements
- add ridge regression in fastglmm and fastlmm
- in variance partitioning analysis, add faster approximation of
  baseline rate for count models

## fastglmm 0.3.6

- add isCountModel() in C++

## fastglmm 0.3.5

- Dec 8, 2025
- Fixed but in C++ code for fastglmm for mu and residuals

## fastglmm 0.3.4

- Nov 24, 2025
- update
  [`varpart()`](http://gabrielhoffman.github.io/fastglmm/reference/varpart.md)
  for NB models
- [`fastglmm()`](http://gabrielhoffman.github.io/fastglmm/reference/fastglmm-pkg.md)
  and [`glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html)
- handle singular models for NB fit
- add Cox-Reid option for NB models
- looser convergence criteria for GLMM

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
  [`fastglmm()`](http://gabrielhoffman.github.io/fastglmm/reference/fastglmm-pkg.md)
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
  [`fastglmm()`](http://gabrielhoffman.github.io/fastglmm/reference/fastglmm-pkg.md)

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
