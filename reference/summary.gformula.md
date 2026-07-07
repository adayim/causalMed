# Summary

Summary method for objects returned by
[`gformula`](https://adayim.github.io/causalMed/reference/gformula.md)
or
[`mediation`](https://adayim.github.io/causalMed/reference/mediation.md).
Shows the full estimation results followed by fitted model coefficient
tables.

## Usage

``` r
# S3 method for class 'gformula'
summary(object, digits = max(3, getOption("digits") - 3), ...)
```

## Arguments

- object:

  Object of class `"gformula"`.

- digits:

  Integer. Number of decimal places used when rounding numeric output.
  Default `max(3, getOption("digits") - 3)`.

- ...:

  Not used.
