# Log Moments of NB given X, Beta, and offset

Log Moments of NB given X, Beta, and offset

## Usage

``` r
log_moments_nb_BlupXBZ(
  BLUP,
  X,
  Beta,
  Z,
  weights,
  offset,
  theta,
  delta,
  method,
  dcmpMethod,
  c = 1,
  p_tail = 1e-04,
  nthreads = 10L
)
```

## Arguments

- BLUP:

  BLUP

- X:

  design matrix

- Beta:

  coefs

- Z:

  random effects design matrix

- weights:

  sample-level weights

- offset:

  offset

- theta:

  overdispersion parameters

- delta:

  ratio of variance components

- method:

  method

- dcmpMethod:

  dcmpMethod

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
library(fastglmm)

data(PsychAD)

# regression formula
form <- PTPRG ~ (1|SubID) + offset(log(libSize))

# NB GLMM on PTPRG expression via PQL
fit <- fastglmm.nb(form, PsychAD)

log_moments_nb_BlupXBZ(
  BLUP = ranef(fit)$SubID,
  X = model.matrix(fit), 
  Beta = as.matrix(coef(fit)), 
  Z = as.matrix(fit$Z),
  weights = fit$prior.weights,
  offset = fit$offset, 
  theta = getTheta(fit),
  delta = fit$delta,
  method = "exact", 
  dcmpMethod = "categorical"
  )
#> Error: unable to find an inherited method for function ‘getTheta’ for signature ‘object = "fastglmm"’
```
