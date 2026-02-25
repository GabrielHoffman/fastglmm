# Object Summaries and Hypothesis Testing

Object summaries and hypothesis testing of fixed effects

## Usage

``` r
# S3 method for class 'fastlmm'
summary(object, ddf = c("satterthwaite", "asymptotic"), ...)
```

## Arguments

- object:

  fitted model of class `fastlmm`

- ddf:

  `"satterthwaite"`: use Satterthwaite approximation to denominator
  degrees of freedom for the Student-t distribution, or `"asymptotic"`
  to use normal distribution as null for the test statistic

- ...:

  other args, not used
