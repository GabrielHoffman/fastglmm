# Negative Binomial Mixed Model

Recent advances in single cell transcriptomics have enabled generation
of large-scale datasets of millions of cells from hundreds or thousands
of samples [(Rood, et al.,
2024)](https://doi.org/10.1038/s41586-024-08338-4). In differential
expression analysis, the goal is to identify genes whose expression is
associated with a variable of interest. Single cell transcriptomics
produces counts for each gene and each observed cell.

To be concrete, consider counts from the gene PTPRG from 249 subjects
measured across 60k microglia cells where an average of 202 cells were
observed per subject [(Hoffman, et al,
2025)](https://doi.org/10.1101/2023.03.17.533005). Consider the dataset
`PsychAD` storing observed cells as rows where column `PTPRG` is the
number of counts for this gene, `libSize` is the total counts for each
cell, `Dx` is Alzheimer’s disease status, `Age` and `Sex` indicate
subject age and sex, and `SubID` indicates which subject each observed
cell is from. Here we use a negative binomial mixed model (NBMM) to test
if the counts of the PTPRG gene are associated with Alzheimer’s disease
status after accounting for the total number of counts for each cell.
Age and sex are included as covariates and subject is a random effect to
account for the repeated measures design. We use the familiar formula
syntax, set `log(libSize)` to be an `offset` term, and use the
[`fastglmm.nb()`](http://gabrielhoffman.github.io/fastglmm/reference/fastglmm.nb.md)
function:

## Fit NBMM

``` r

library(fastglmm)

data(PsychAD)

form <- PTPRG ~ offset(log(libSize)) + Dx + Age + Sex + (1|SubID)
fit <- fastglmm.nb(form, PsychAD)

fit
```

    ## 
    ## Call:
    ## fastglmm(formula = formula, data = data, family = negative.binomial(NA), 
    ##     weights = weights, maxit = maxit, tol = tol, tol.eta = tol.eta, 
    ##     doCoxReid = doCoxReid, nthreads = nthreads)
    ## 
    ## Coefficients:
    ## (Intercept)         DxAD          Age      SexMale  
    ##   -8.577681     1.102755    -0.007262    -0.103277

## Hypothesis testing

``` r

summary(fit)
```

    ## Generalized linear mixed model fit by PQL ['fastglmm']
    ##  Family: Negative Binomial(0.2443)  ( log )
    ##  Formula: PTPRG ~ offset(log(libSize)) + Dx + Age + Sex + (1 | SubID)
    ## 
    ## Coefficients:
    ##               Estimate Std. Error    df t value Pr(>|t|)    
    ## (Intercept)  -8.577681   0.311172 149.9 -27.566   <2e-16 ***
    ## DxAD          1.102755   0.073541 145.6  14.995   <2e-16 ***
    ## Age          -0.007262   0.003775 148.2  -1.923   0.0563 .  
    ## SexMale      -0.103277   0.075756 145.7  -1.363   0.1749    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual df: 60227.2 
    ## 
    ## Variance components:
    ##   sigSq_g: 0.4803
    ##   sigSq_e: 8.17
    ##   delta:   17.01

The PTPRG gene has increased expression in subjects with Alzheimer’s
disease compared with controls with log fold change 1.10277 and a
p-value \< 2\times 10^{-16}. The count data shows high overdispersion
with \theta = 0.2443, where \theta = \infty is equivalent to a Poisson
distribution with no overdispersion.

Note that while 71% of the 60k observations are zero, this is due to a
low count rate rather than zero inflation. This is supported by
comparison with the zero inflation component using
[`glmmTMB::glmmTMB()`](https://rdrr.io/pkg/glmmTMB/man/glmmTMB.html)
(p-value of zero-inflation intercept: 0.961), and consistent with other
work on single cell transcriptomics ([Jiang, et al.,
2022](https://doi.org/10.1186/s13059-022-02601-5), [Sarkar and Stephens,
2021](https://doi.org/10.1038/s41588-021-00873-4), [Svensson,
2020](https://doi.org/10.1038/s41587-019-0379-5)).

## Variance partitioning analysis

``` r

varpart(fit)
```

    ##           Dx          Age          Sex        SubID   CountNoise    Residuals 
    ## 0.0133097371 0.0002197195 0.0001160483 0.0210313340 0.3062036641 0.6591194970

Variance partitioning analysis quantifies the contribution of each
variable to the observed variance in PTPRG expression. Alzheimer’s
disease status explains 1.3% of observed expression variation, variation
across subjects explains 2.1%, with sex and age making a smaller
contribution. Measurement error due to sampling a finite number of
counts (i.e. `CountNoise`) explains 30.6% of the variance, and the
remaining 65.9% of the variance is explained by the residuals.

##### Session info

``` r

sessionInfo()
```

    ## R version 4.5.1 (2025-06-13)
    ## Platform: aarch64-apple-darwin23.6.0
    ## Running under: macOS Sonoma 14.7.1
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /opt/homebrew/Cellar/openblas/0.3.31_1/lib/libopenblasp-r0.3.31.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/New_York
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] fastglmm_0.4.1 nlme_3.1-168  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] sass_0.4.10       generics_0.1.4    lattice_0.22-9    lme4_2.0-1       
    ##  [5] digest_0.6.39     magrittr_2.0.4    evaluate_1.0.5    grid_4.5.1       
    ##  [9] fastmap_1.2.0     jsonlite_2.0.0    Matrix_1.7-4      Formula_1.2-5    
    ## [13] codetools_0.2-20  textshaping_1.0.5 jquerylib_0.1.4   reformulas_0.4.4 
    ## [17] abind_1.4-8       Rdpack_2.6.6      cli_3.6.5         rlang_1.1.7      
    ## [21] rbibutils_2.4.1   splines_4.5.1     cachem_1.1.0      yaml_2.3.12      
    ## [25] otel_0.2.0        tools_4.5.1       nloptr_2.2.1      minqa_1.2.8      
    ## [29] dplyr_1.2.0       boot_1.3-32       vctrs_0.7.1       R6_2.6.1         
    ## [33] matrixStats_1.5.0 lifecycle_1.0.5   fs_1.6.7          car_3.1-5        
    ## [37] htmlwidgets_1.6.4 MASS_7.3-65       ragg_1.5.1        pkgconfig_2.0.3  
    ## [41] desc_1.4.3        pkgdown_2.2.0     pillar_1.11.1     bslib_0.10.0     
    ## [45] glue_1.8.0        Rcpp_1.1.1        systemfonts_1.3.2 xfun_0.56        
    ## [49] tibble_3.3.1      tidyselect_1.2.1  knitr_1.51        htmltools_0.5.9  
    ## [53] rmarkdown_2.30    carData_3.0-6     compiler_4.5.1
