# Create sparse indicator matrix

Create sparse indicator matrix from factor

## Usage

``` r
preprocess_indicator(x)
```

## Arguments

- x:

  a `factor`

## Value

sparse indicator matrix with levels as columns

## Examples

``` r
data(PsychAD)

# extract a few subjects
idx <- c(1, 4, 100, 444)
dcmp <- preprocess_indicator(PsychAD$SubID[idx])
```
