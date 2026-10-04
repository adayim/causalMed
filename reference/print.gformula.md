# Print the results of gformula() or mediation()

Prints the estimated mean outcome (or risk) under each intervention, the
contrasts or mediation decomposition, a legend of the interventions, the
analysis setup, and the observed nonparametric benchmark.

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
