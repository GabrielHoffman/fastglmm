# Gene expression of PTPRG in microglia in 60k cells

A dataset containing gene expression of PTPRG in 60k microglia cells
from the PsychAD Consortium from 299 subjects. Of these, 149 subjects
have Alzheimer's disease, and 150 are neurotypical controls.

## Usage

``` r
data(PsychAD)
```

## Format

A data frame with 60481 rows (i.e. cells) and 6 variables:

- SubID:

  subject identifier

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

## Value

none

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
#>  Family: Negative Binomial(0.2413)  ( log )
#>  Formula: PTPRG ~ Dx + Age + Sex + offset(log(libSize)) + (1 | SubID)
#> 
#> Coefficients:
#>               Estimate Std. Error    df t value Pr(>|t|)    
#> (Intercept)  -8.579598   0.314805 149.9 -27.254   <2e-16 ***
#> DxAD          1.103112   0.074407 145.6  14.825   <2e-16 ***
#> Age          -0.007246   0.003820 148.2  -1.897   0.0598 .  
#> SexMale      -0.102996   0.076649 145.7  -1.344   0.1811    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual df: 60226.7 
#> 
#> Variance components:
#>   sigSq_g: 0.4883
#>   sigSq_e: 8.176
#>   delta:   16.74

# Variance partitioning analysis
varpart(fit)
#>           Dx          Age          Sex        SubID   CountNoise    Residuals 
#> 0.0375459290 0.0006167199 0.0003253761 0.0602769180 0.3940102045 0.5072248524 
```
