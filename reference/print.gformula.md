# Print

Print method for objects returned by
[`gformula`](https://adayim.github.io/causalMed/reference/gformula.md)
or
[`mediation`](https://adayim.github.io/causalMed/reference/mediation.md).

## Usage

``` r
# S3 method for class 'gformula'
print(x, models = FALSE, digits = max(3, getOption("digits") - 3), ...)
```

## Arguments

- x:

  Object of class `"gformula"`.

- models:

  Logical. If `TRUE`, print fitted model details (call and coefficients,
  or full summary when `return_fitted = TRUE`). Default `FALSE`.

- digits:

  Integer. Number of decimal places used when rounding numeric output.
  Default `max(3, getOption("digits") - 3)`.

- ...:

  Not used.
