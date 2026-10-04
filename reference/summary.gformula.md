# Summarise the results of gformula() or mediation()

Prints the same results as
[`print.gformula`](https://adayim.github.io/causalMed/reference/print.gformula.md),
followed by the coefficient table of every fitted model.

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
