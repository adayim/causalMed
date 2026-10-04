# Error Catching

This function is used to check for errors in the
[`gformula`](https://adayim.github.io/causalMed/reference/gformula.md).

## Usage

``` r
check_error(data, id_var, base_vars, exposure, time_var, models)
```

## Arguments

- data:

  A `data.frame` in long format.

- id_var:

  Character. Name of the subject identifier.

- base_vars:

  Character vector of time-fixed baseline covariates (may be empty),
  with no missing values. Only these and `id_var` are carried into the
  simulated cohort; every other variable a model uses, including lags,
  must be created by the recode hooks.

- exposure:

  Character. Name of the exposure to intervene on.

- time_var:

  Character. Name of the numeric time variable. Each distinct value is
  one simulated step, in ascending order.

- models:

  List of
  [`spec_model`](https://adayim.github.io/causalMed/reference/spec_model.md)
  objects, in the order the variables are generated.

## Value

No value is returned.
