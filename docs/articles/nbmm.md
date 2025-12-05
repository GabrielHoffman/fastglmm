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
    ##     nthreads = nthreads)
    ## 
    ## Coefficients:
    ## (Intercept)         DxAD          Age      SexMale  
    ##   -8.577764     1.102771    -0.007261    -0.103264

## Hypothesis testing

``` r

summary(fit)
```

    ## Generlized linear mixed model fit by PQL ['fastglmm']
    ##  Family: Negative Binomial(0.2443)  ( log )
    ##  Formula: PTPRG ~ offset(log(libSize)) + Dx + Age + Sex + (1 | SubID)
    ## 
    ## Coefficients:
    ##              Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -8.577764   0.369236 -23.231   <2e-16 ***
    ## DxAD         1.102771   0.087264  12.637   <2e-16 ***
    ## Age         -0.007261   0.004480  -1.621    0.105    
    ## SexMale     -0.103264   0.089893  -1.149    0.251    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual df: 60227.2 
    ## 
    ## Variance components:
    ##   sigSq_g: 0.481
    ##   sigSq_e: 8.178
    ##   delta:   17

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

    ##           Dx          Age          Sex        SubID    Residuals 
    ## 0.0133011914 0.0002195276 0.0001159415 0.0210472556 0.9653160840

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
    ## BLAS/LAPACK: /opt/homebrew/Cellar/openblas/0.3.30/lib/libopenblasp-r0.3.30.dylib;  LAPACK version 3.12.0
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
    ## [1] fastglmm_0.3.3 lme4_1.1-37    Matrix_1.7-4  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] sass_0.4.10         generics_0.1.4      lattice_0.22-7     
    ##  [4] digest_0.6.37       magrittr_2.0.4      evaluate_1.0.5     
    ##  [7] grid_4.5.1          fastmap_1.2.0       jsonlite_2.0.0     
    ## [10] Formula_1.2-5       numDeriv_2016.8-1.1 textshaping_1.0.4  
    ## [13] jquerylib_0.1.4     abind_1.4-8         reformulas_0.4.2   
    ## [16] Rdpack_2.6.4        cli_3.6.5           rlang_1.1.6        
    ## [19] rbibutils_2.4       splines_4.5.1       cachem_1.1.0       
    ## [22] yaml_2.3.10         tools_4.5.1         nloptr_2.2.1       
    ## [25] minqa_1.2.8         dplyr_1.1.4         boot_1.3-32        
    ## [28] vctrs_0.6.5         R6_2.6.1            matrixStats_1.5.0  
    ## [31] lifecycle_1.0.4     fs_1.6.6            car_3.1-3          
    ## [34] htmlwidgets_1.6.4   MASS_7.3-65         ragg_1.5.0         
    ## [37] pkgconfig_2.0.3     desc_1.4.3          pkgdown_2.2.0      
    ## [40] bslib_0.9.0         pillar_1.11.1       glue_1.8.0         
    ## [43] Rcpp_1.1.0          systemfonts_1.3.1   xfun_0.54          
    ## [46] tibble_3.3.0        tidyselect_1.2.1    knitr_1.50         
    ## [49] htmltools_0.5.8.1   nlme_3.1-168        rmarkdown_2.30     
    ## [52] carData_3.0-5       compiler_4.5.1
