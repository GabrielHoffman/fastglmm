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
