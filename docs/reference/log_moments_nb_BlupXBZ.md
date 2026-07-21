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
