# Variance Partitioning Analysis

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`fastglmm`](https://github.com/GabrielHoffman/fastlmm)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`mvtnorm`](http://mvtnorm.R-forge.R-project.org)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`MASS`](http://www.stats.ox.ac.uk/pub/MASS4/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`glmmTMB`](https://github.com/glmmTMB/glmmTMB)`)`

## Simulation for Gaussian regression

\
`sim_data`` ``<-`` ``function``(``n``, ``beta``, ``sigSq_g``, ``rho``, ``n_random``, ``intercept``, ``link``, ``signal_frac`` ``=`` ``NULL``, ``theta``=``NULL``)``{`\
\
`  ``# Create fixed covariates with given correlation rho`\
`  ``Sigma`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``1``,``rho``,``rho``,``1``)``, ``2``)`\
`  ``X`` ``<-`` `[`scale`](https://rdrr.io/r/base/scale.html)`(`[`rmvnorm`](https://rdrr.io/pkg/mvtnorm/man/Mvnorm.html)`(``n``, `[`c`](https://rdrr.io/r/base/c.html)`(``0``,``0``)``, ``Sigma``)``)`\
\
`  ``# random effecdt with n_random levels`\
`  ``z`` ``<-`` `[`factor`](https://rdrr.io/r/base/factor.html)`(`[`sample`](https://rdrr.io/r/base/sample.html)`(`[`seq`](https://rdrr.io/r/base/seq.html)`(``n_random``)``, ``n``, replace``=``TRUE``)``)`\
\
`  ``# Sparse design matrix for random effects`\
`  ``df`` ``<-`` `[`data.frame`](https://rdrr.io/r/base/data.frame.html)`(``X``, ``z``)`\
`  ``Z`` ``<-`` `[`sparse.model.matrix`](https://rdrr.io/pkg/Matrix/man/sparse.model.matrix.html)`(``~`` ``0`` ``+`` ``z``, ``df``)`\
\
`  ``# Simulate random effect`\
`  ``alpha`` ``<-`` `[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(`[`ncol`](https://rdrr.io/r/base/nrow.html)`(``Z``)``, ``0``, sd``=`[`sqrt`](https://rdrr.io/r/base/MathFun.html)`(``sigSq_g``)``)`\
\
`  ``# linear predictor based on interceptm, fixed & random effects`\
`  ``eta`` ``<-`` ``intercept`` ``+`` `[`as.matrix`](https://rdrr.io/r/base/matrix.html)`(``X`` `[`%*%`](https://rdrr.io/r/base/matmult.html)` ``beta`` ``+`` ``Z`` `[`%*%`](https://rdrr.io/r/base/matmult.html)` ``alpha``)`\
\
`  ``# noise in Gaussian model`\
`  ``sigSq_e`` ``<-`` ``(``1``-``signal_frac``)``/``signal_frac`` ``*`` `[`var`](https://rdrr.io/r/stats/cor.html)`(``eta``)``[``1``]`\
\
`  ``# Simulate response`\
`  ``df``$``y`` ``<-`` `[`switch`](https://rdrr.io/r/base/switch.html)`(``link``, `\
`    ``"gaussian"`` ``=`` ``eta`` ``+`` `[`rnorm`](https://rdrr.io/r/stats/Normal.html)`(``n``, ``0``, sd ``=`` `[`sqrt`](https://rdrr.io/r/base/MathFun.html)`(``sigSq_e``)``)``,`\
`    ``"poisson"`` ``=`` `[`rpois`](https://rdrr.io/r/stats/Poisson.html)`(``n``, `[`exp`](https://rdrr.io/r/base/Log.html)`(``eta``)``)``,`\
`    ``"logistic"`` ``=`` `[`rbinom`](https://rdrr.io/r/stats/Binomial.html)`(``n``, ``1``, prob``=`[`plogis`](https://rdrr.io/r/stats/Logistic.html)`(``eta``)``)``,`\
`    ``"nb"`` ``=`` `[`rnegbin`](https://rdrr.io/pkg/MASS/man/rnegbin.html)`(``n``, `[`exp`](https://rdrr.io/r/base/Log.html)`(``eta``)``, theta ``=`` ``theta``)``)`\
\
`  ``# Unexplained variance for each model`\
`  ``sigSq_noise`` ``<-`` `[`switch`](https://rdrr.io/r/base/switch.html)`(``link``, `\
`    ``"gaussian"`` ``=`` ``sigSq_e``,`\
`    ``"poisson"`` ``=``  `[`trigamma`](https://rdrr.io/r/base/Special.html)`(`[`exp`](https://rdrr.io/r/base/Log.html)`(``intercept``)``)``,`\
`    ``"logistic"`` ``=`` ``pi``^``2`` ``/`` ``3``,`\
`    ``"nb"`` ``=`` `[`trigamma`](https://rdrr.io/r/base/Special.html)`(``1``/``(``1``/`[`exp`](https://rdrr.io/r/base/Log.html)`(``intercept``)`` ``+`` ``1``/``theta``)``)``)`\
\
`  ``# variance components`\
`  ``# varTotal = `\
`  ``# vc = c(beta^2, sigSq_g, sigSq_noise)`\
\
`  ``signal_frac`` ``=`` `[`var`](https://rdrr.io/r/stats/cor.html)`(``eta``)``[``1``]`` ``/`` ``(`[`var`](https://rdrr.io/r/stats/cor.html)`(``eta``)``[``1``]`` ``+`` ``sigSq_noise``)`\
\
`  ``a`` ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``beta``^``2``, ``sigSq_g``)`\
`  ``vc`` ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``signal_frac`` ``*`` ``a``/`[`sum`](https://rdrr.io/r/base/sum.html)`(``a``)``, ``1`` ``-`` ``signal_frac``)`\
\
`  `[`attr`](https://rdrr.io/r/base/attr.html)`(``df``, ``"varFrac"``)`` ``=`` ``vc`` ``/`` `[`sum`](https://rdrr.io/r/base/sum.html)`(``vc``)`\
`  ``df`\
`}`\
\
\
`df`` ``<-`` ``sim_data``(``n ``=`` ``1e6``, `\
`              beta ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``1``,``2``)``, `\
`              sigSq_g ``=`` ``2``,`\
`              rho ``=`` ``.7``, `\
`              n_random ``=`` ``1000``,`\
`              intercept ``=`` ``3``,`\
`              ``"gaussian"``,`\
`              signal_frac ``=`` ``.8``)`\
\
\
`fit`` ``<-`` `[`fastlmm`](http://gabrielhoffman.github.io/fastglmm/reference/fastlmm.md)`(``y`` ``~`` ``X1`` ``+`` ``X2`` ``+`` ``(``1``|``z``)``, ``df``)`\
\
[`varpart`](http://gabrielhoffman.github.io/fastglmm/reference/varpart.md)`(``fit``)`

    ##        X1        X2         z Residuals 
    ## 0.1159611 0.4638733 0.2202469 0.1999187

\
[`attr`](https://rdrr.io/r/base/attr.html)`(``df``, ``"varFrac"``)`

    ## [1] 0.1142857 0.4571429 0.2285714 0.2000000

## Fixed intercept for Z, now fix

## true variance fractions for non-normal

\
`df`` ``<-`` ``sim_data``(``n ``=`` ``1e6``, `\
`              beta ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``1``,``2``)``, `\
`              sigSq_g ``=`` ``1``,`\
`              rho ``=`` ``.7``, `\
`              n_random ``=`` ``500``,`\
`              intercept ``=`` ``-``3``,`\
`              ``"poisson"``)`\
\
[`hist`](https://rdrr.io/r/graphics/hist.html)`(``df``$``y``)`

![](variancePartitioning_files/figure-html/poisson-1.png)

\
`fit`` ``<-`` `[`fastglmm`](http://gabrielhoffman.github.io/fastglmm/reference/fastglmm.md)`(``y`` ``~`` ``X1`` ``+`` ``X2`` ``+`` ``(``1``|``z``)``, ``df``, family``=`[`poisson`](https://rdrr.io/r/stats/family.html)`(``)``)`\
\
[`varpart`](http://gabrielhoffman.github.io/fastglmm/reference/varpart.md)`(``fit``)`

    ##        X1        X2         z Residuals 
    ## 0.1601416 0.6402481 0.1725216 0.0270887

\
[`attr`](https://rdrr.io/r/base/attr.html)`(``df``, ``"varFrac"``)`

    ## [1] 0.00357594 0.01430376 0.00357594 0.97854436

\
`it`` ``=`` `[`glmmTMB`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html)`(``y`` ``~`` ``X1`` ``+`` ``X2`` ``+`` ``(``1``|``z``)``, ``df``, family``=`[`poisson`](https://rdrr.io/r/stats/family.html)`(``)``)`\
\
\
`# No variance component`\
`df`` ``<-`` ``sim_data``(``n ``=`` ``1e6``, `\
`              beta ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``1``,``.2``)``, `\
`              sigSq_g ``=`` ``1e-5``,`\
`              rho ``=`` ``.7``, `\
`              n_random ``=`` ``500``,`\
`              intercept ``=`` ``2``,`\
`              ``"poisson"``)`\
\
[`hist`](https://rdrr.io/r/graphics/hist.html)`(``df``$``y``)`

![](variancePartitioning_files/figure-html/poisson-2.png)

\
`fit`` ``<-`` `[`fastglmm`](http://gabrielhoffman.github.io/fastglmm/reference/fastglmm.md)`(``y`` ``~`` ``X1`` ``+`` ``X2`` ``+`` ``(``1``|``z``)``, ``df``, family``=`[`poisson`](https://rdrr.io/r/stats/family.html)`(``)``)`\
\
[`varpart`](http://gabrielhoffman.github.io/fastglmm/reference/varpart.md)`(``fit``)`

    ##           X1           X2            z    Residuals 
    ## 9.115835e-01 3.629820e-02 1.161418e-05 5.210670e-02

\
[`attr`](https://rdrr.io/r/base/attr.html)`(``df``, ``"varFrac"``)`

    ## [1] 8.664114e-01 3.465646e-02 8.664114e-06 9.892343e-02

\
[`median`](https://rdrr.io/r/stats/median.html)`(``df``$``y``)`

    ## [1] 7

\
`fit1`` ``=`` `[`glm`](https://rdrr.io/r/stats/glm.html)`(``y`` ``~`` ``X1`` ``+`` ``X2``, ``df``, family``=`[`poisson`](https://rdrr.io/r/stats/family.html)`(``)``)`\
[`varpart`](http://gabrielhoffman.github.io/fastglmm/reference/varpart.md)`(``fit1``)`

    ##         X1         X2  Residuals 
    ## 0.91159427 0.03629590 0.05210983

\
`# HOW DOES RSQ COMPARE????`\
`rsq``::`[`rsq`](https://rdrr.io/pkg/rsq/man/rsq.html)`(``fit1``, type``=``"v"``)`

    ## [1] 0.9748794

\
`w`` ``=`` ``1``/``(``df``$``y`` ``+`` ``.1``)`\
`fit1`` ``=`` `[`lm`](https://rdrr.io/r/stats/lm.html)`(`[`log`](https://rdrr.io/r/base/Log.html)`(``y`` ``+`` ``.1``)`` ``~`` ``X1`` ``+`` ``X2``, ``df``, weights``=``w`` ``/`` `[`mean`](https://rdrr.io/r/base/mean.html)`(``w``)``)`\
[`varpart`](http://gabrielhoffman.github.io/fastglmm/reference/varpart.md)`(``fit1``)`

    ##         X1         X2  Residuals 
    ## 0.51228583 0.01960208 0.46811210

## Application to PsychAD dataset

\
[`data`](https://rdrr.io/r/utils/data.html)`(``PsychAD``)`\
\
`form`` ``=`` ``PTPRG`` ``~`` ``Dx`` ``+`` ``Age`` ``+`` ``Sex`` ``+`` `[`offset`](https://rdrr.io/r/stats/offset.html)`(`[`log`](https://rdrr.io/r/base/Log.html)`(``libSize``)``)`` ``+`` ``(``1``|``SubID``)`\
\
`fit`` ``<-`` `[`fastglmm.nb`](http://gabrielhoffman.github.io/fastglmm/reference/fastglmm.nb.md)`(``form``, ``PsychAD``)`\
\
[`varpart`](http://gabrielhoffman.github.io/fastglmm/reference/varpart.md)`(``fit``)`

    ##           Dx          Age          Sex        SubID    Residuals 
    ## 0.0133011913 0.0002195276 0.0001159415 0.0210472556 0.9653160840
