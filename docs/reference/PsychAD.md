# Gene expession of PTPRG in microglia in 60k cells

A dataset containing gene expression of PTPRG in 60k microglia cells
from the PsychAD Consortium from 299 subjects. Of these, 149 subjects
have Alzheimer's disease, and 150 are neurotypical controls.

## Usage

``` r
PsychAD
```

## Format

A data frame with 60481 rows (i.e. cells) and 6 variables:

- SubID:

  subject identifer

- Sex:

  Sex

- Dx:

  Diagnosis: AD or Control

- Age:

  subject age

- PTPRG:

  number of RNA-seq counts for this gene

- libSize:

  total number of RNA-seq reads for this cell

## Source

[doi:10.1101/2023.03.17.533005](https://doi.org/10.1101/2023.03.17.533005)
, <https://diseaseneurogenomics.github.io/dreamlet/index.html>

## Examples

``` r
data(PsychAD)

# model formula
form <- PTPRG ~ Dx + Age + Sex + offset(log(libSize)) + (1|SubID)

# fit negative binomial mixed model
fam <- negative.binomial(NA)
fit <- fastglmm(form, PsychAD, family=fam)

summary(fit)
#> Generalized linear mixed model fit by PQL ['fastglmm']
#>  Family: Negative Binomial(0.2443)  ( log )
#>  Formula: PTPRG ~ Dx + Age + Sex + offset(log(libSize)) + (1 | SubID)
#> 
#> Coefficients:
#>               Estimate Std. Error   ddf t value Pr(>|t|)    
#> (Intercept)  -8.577681   0.311172 149.9 -27.566   <2e-16 ***
#> DxAD          1.102755   0.073541 145.6  14.995   <2e-16 ***
#> Age          -0.007262   0.003775 148.2  -1.923   0.0563 .  
#> SexMale      -0.103277   0.075756 145.7  -1.363   0.1749    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual df: 60227.2 
#> 
#> Variance components:
#>   sigSq_g: 0.4803
#>   sigSq_e: 8.17
#>   delta:   17.01

# Variance partitioning analysis
varpart(fit)
#>           Dx          Age          Sex        SubID   CountNoise    Residuals 
#> 0.0133097371 0.0002197195 0.0001160483 0.0210313340 0.3062036641 0.6591194970 
```
