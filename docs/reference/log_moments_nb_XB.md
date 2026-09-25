# Log Moments of NB given X, Beta, and offset

Log Moments of NB given X, Beta, and offset

## Usage

``` r
log_moments_nb_XB(
  X,
  Beta,
  offset,
  theta,
  method,
  c = 1,
  p_tail = 1e-04,
  nthreads = 10L
)
```

## Arguments

- X:

  design matrix

- Beta:

  coefs

- offset:

  offset

- theta:

  overdispersion parameters

- method:

  method

- c:

  pseudocount

- p_tail:

  probability cutoff

- nthreads:

  number of threads

## Value

variance of the signal (var.signal), variance of the noise (var.noise),
(alpha) fraction of Poisson noise

## Examples

``` r
library(MASS)
data(PsychAD)

# regression formula
form <- PTPRG ~ offset(log(libSize))

if (FALSE) { # \dontrun{
# NB GLM on PTPRG expression 
fit <- glm.nb(form, PsychAD)

log_moments_nb_XB(
  model.matrix(fit), 
  as.matrix(coef(fit)), 
  fit$offset, 
  getNBTheta(fit), 
  method="exact")
} # }
```
