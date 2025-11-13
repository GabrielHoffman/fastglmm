# Causal Mediation Analysis for fastglmm

Estimate various quantities for causal mediation analysis, including
average causal mediation effects (indirect effect), average direct
effects, proportions mediated, and total effect. Here, adapted from
[`mediation::mediate()`](https://rdrr.io/pkg/mediation/man/mediate.html)
to handle models fit with
[`fastlmm()`](http://gabrielhoffman.github.io/fastglmm/reference/fastlmm.md)
and
[`fastglmm()`](http://gabrielhoffman.github.io/fastglmm/reference/fastglmm.md).

## Usage

``` r
mediate(model.m, model.y, sims = 1000, treat, mediator)
```

## Arguments

- model.m:

  a fitted model object for mediator.

- model.y:

  a fitted model object for outcome.

- sims:

  number of Monte Carlo draws for quasi-Bayesian approximation

- treat:

  a character string indicating the name of the treatment variable used
  in the models. The treatment can be either binary (integer or a
  two-valued factor) or continuous (numeric).

- mediator:

  a character string indicating the name of the mediator variable used
  in the models.

## See also

[`mediation::mediate()`](https://rdrr.io/pkg/mediation/man/mediate.html)

## Examples

``` r
# Model
# x = 2*z + noise
# y = 0*x - 1*z + noise

n <- 1000
data <- data.frame(z = rnorm(n))
data$x <- with(data, z * 2 + rnorm(n,0,3))
data$y <- with(data, x * 0 + -1*z + rnorm(n,0,3))
data$letter <- sample(LETTERS[1:2], n, replace=TRUE)

fit1 <- glm(x ~ z, data=data)
fit2 <- glm(y ~ x + z, data=data)

# Estimation via quasi-Bayesian approximation
contcont <- mediation::mediate(fit1, fit2, sims=1000, treat="z", mediator="x")
summary(contcont)
#> 
#> Causal Mediation Analysis 
#> 
#> Quasi-Bayesian Confidence Intervals
#> 
#>                 Estimate 95% CI Lower 95% CI Upper p-value    
#> ACME           -0.038269    -0.172673     0.101577    0.57    
#> ADE            -0.783406    -1.014122    -0.558230  <2e-16 ***
#> Total Effect   -0.821674    -0.999140    -0.645053  <2e-16 ***
#> Prop. Mediated  0.046671    -0.122533     0.213947    0.57    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Sample Size Used: 1000 
#> 
#> 
#> Simulations: 1000 
#> 

# Include random effect with fastlmm
fit1m <- fastlmm(x ~ z + (1|letter), data=data)
fit2m <- fastlmm(y ~ x + z + (1|letter), data=data)

mediate(fit1m, fit2m, sims=100, treat="z", mediator="x")
#> 
#> Causal Mediation Analysis: fastglmm
#> 
#> Quasi-Bayesian Confidence Intervals
#> 
#>              Estimate       se   CI.low  CI.high  p.value    
#> ACME         -0.04219  0.06938 -0.17729  0.07721    0.543    
#> ADE          -0.79635  0.12413 -1.01788 -0.60918 1.40e-10 ***
#> TotalEffect  -0.83854  0.10219 -1.00972 -0.68639 2.29e-16 ***
#> PropMediated  0.05123  0.08436 -0.09184  0.20357    0.544    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
```
