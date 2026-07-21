# Package index

## Core functions

- [`fastlmm()`](http://gabrielhoffman.github.io/fastglmm/reference/fastlmm.md)
  : Fast Linear Mixed Model with 1 Random Effect
- [`fastglmm()`](http://gabrielhoffman.github.io/fastglmm/reference/fastglmm.md)
  : Fit Generalized Linear Mixed Model via PQL
- [`fastglmm.nb()`](http://gabrielhoffman.github.io/fastglmm/reference/fastglmm.nb.md)
  : Fit Negative Binomial Mixed Model via PQL

## Example Data

- [`PsychAD`](http://gabrielhoffman.github.io/fastglmm/reference/PsychAD.md)
  : Gene expession of PTPRG in microglia in 60k cells

## Generics

- [`anova(`*`<fastlmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/anova.md)
  : ANOVA Tables
- [`coef(`*`<fastlmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/coef.fastlmm.md)
  : Extract Model Coefficients
- [`cooks.distance(`*`<fastlmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/cooks.distance.fastlmm.md)
  : Cook's Distance Metric
- [`deviance(`*`<fastlmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/deviance.fastlmm.md)
  : Model Deviance
- [`ddf()`](http://gabrielhoffman.github.io/fastglmm/reference/ddf.md) :
  Denominator Degrees of Freedom
- [`df.residual(`*`<fastlmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/df.residual.fastlmm.md)
  : Residual Degrees-of-Freedom
- [`edf()`](http://gabrielhoffman.github.io/fastglmm/reference/edf.md) :
  Effective Degrees-of-Freedom of Model Fit
- [`extractAIC(`*`<fastlmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/extractAIC.fastlmm.md)
  : Extract AIC from a Fitted Model
- [`family(`*`<fastlmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/family.md)
  [`family(`*`<fastglmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/family.md)
  [`family(`*`<glmmPQL>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/family.md)
  : Extract Model Family
- [`fitted(`*`<fastlmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/fitted.md)
  [`fitted(`*`<fastglmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/fitted.md)
  : Extract Model Fitted Values
- [`fixef(`*`<fastlmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/fixef.fastlmm.md)
  : Extract Fixed Effects
- [`hatvalues(`*`<fastlmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/hatvalues.fastlmm.md)
  : Return Diagonals of Hat Matrix
- [`linearHypothesis()`](http://gabrielhoffman.github.io/fastglmm/reference/linearHypothesis.md)
  : Test Linear Hypothesis
- [`logLik(`*`<fastlmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/logLik.fastlmm.md)
  : Extract Log-Likelihood
- [`model.frame(`*`<fastlmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/model.frame.fastlmm.md)
  : Extracting the Model Frame from a Fit
- [`model.matrix(`*`<fastlmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/model.matrix.fastlmm.md)
  : Construct Design Matrices
- [`nobs(`*`<fastlmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/nobs.fastlmm.md)
  : Extract the Number of Observations from a Fit
- [`plot(`*`<fastlmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/plot.fastlmm.md)
  : Diagnostic Plots for Model Fits
- [`predict(`*`<fastlmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/predict.md)
  [`predict(`*`<fastglmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/predict.md)
  : Model Predictions
- [`ranef(`*`<fastlmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/ranef.md)
  : Extract the modes of the random effect
- [`residuals(`*`<fastlmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/residuals.fastlmm.md)
  : Extract Model Residuals
- [`rstudent(`*`<fastlmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/rstudent.fastlmm.md)
  : Extract Studentized Residuals
- [`sigma(`*`<fastlmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/sigma.md)
  [`sigma(`*`<fastglmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/sigma.md)
  : Extract Residual Standard Deviation
- [`simulate(`*`<fastlmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/simulate.md)
  : Simulate Responses
- [`summary(`*`<fastlmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/summary.fastlmm.md)
  : Object Summaries and Hypothesis Testing
- [`terms(`*`<fastlmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/terms.fastlmm.md)
  : Model Terms
- [`varpart()`](http://gabrielhoffman.github.io/fastglmm/reference/varpart.md)
  : Variance Partitioning Analysis
- [`vcov(`*`<fastlmm>`*`)`](http://gabrielhoffman.github.io/fastglmm/reference/vcov.md)
  : Calculate Variance-Covariance Matrix for a Fitted Model Object

## Power analysis

- [`meff()`](http://gabrielhoffman.github.io/fastglmm/reference/meff.md)
  : Effective Number of Independent Measurements per Subject
- [`plotSaturation()`](http://gabrielhoffman.github.io/fastglmm/reference/plotSaturation.md)
  : Plot Saturation Curve for NBMM
- [`powerNBMM()`](http://gabrielhoffman.github.io/fastglmm/reference/powerNBMM.md)
  : Power for Negative Binomial Mixed Model

## Other

- [`as.fastlmm()`](http://gabrielhoffman.github.io/fastglmm/reference/as.fastlmm.md)
  : Convert list to fastlmm class
- [`getFamilyString()`](http://gabrielhoffman.github.io/fastglmm/reference/getFamilyString.md)
  : Convert GLM family to string
- [`indicator_decomp()`](http://gabrielhoffman.github.io/fastglmm/reference/indicator_decomp.md)
  : Spectral decomposition of factor indicator matrix
- [`preprocess_indicator()`](http://gabrielhoffman.github.io/fastglmm/reference/preprocess_indicator.md)
  : Create sparse indicator matrix
- [`fastlmm.fit()`](http://gabrielhoffman.github.io/fastglmm/reference/fastlmm.fit.md)
  : Fitter Function for Linear Mixed Model
- [`nb_theta()`](http://gabrielhoffman.github.io/fastglmm/reference/nb_theta.md)
  : Estimate theta of the Negative Binomial
- [`refitModel()`](http://gabrielhoffman.github.io/fastglmm/reference/refitModel.md)
  : Refit fastlmm model with new delta
- [`dispersion()`](http://gabrielhoffman.github.io/fastglmm/reference/dispersion.md)
  : Overdispersion parameter phi for quasi-likelihood
- [`varianceTerms()`](http://gabrielhoffman.github.io/fastglmm/reference/varianceTerms.md)
  : Array of variance component estimates
- [`negative.binomial()`](http://gabrielhoffman.github.io/fastglmm/reference/negative.binomial.md)
  : Family function for Negative Binomial GLMs

## Pure R Implementations (for intuition only)

- [`fastlmm_R()`](http://gabrielhoffman.github.io/fastglmm/reference/fastlmm_R.md)
  : Fit linear mixed model using SVD of covariance
- [`fastglmm_R()`](http://gabrielhoffman.github.io/fastglmm/reference/fastglmm_R.md)
  : Fit generalized linear mixed model via PQL
- [`fastglmm_R.nb()`](http://gabrielhoffman.github.io/fastglmm/reference/fastglmm_R.nb.md)
  : Fit negative binomial mixed model via PQL
