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

<https://doi.org/10.1101/2023.03.17.533005>,
<https://diseaseneurogenomics.github.io/dreamlet/index.html>

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
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -8.577763   0.369235 -23.231   <2e-16 ***
#> DxAD         1.102771   0.087264  12.637   <2e-16 ***
#> Age         -0.007261   0.004480  -1.621    0.105    
#> SexMale     -0.103264   0.089892  -1.149    0.251    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Residual df: 60227.2 
#> 
#> Variance components:
#>   sigSq_g: 0.481
#>   sigSq_e: 8.178
#>   delta:   17

# Variance partitioning analysis
varpart(fit)
#>           Dx          Age          Sex        SubID    Residuals 
#> 0.0133012726 0.0002195293 0.0001159424 0.0210471707 0.9653160849 
```
