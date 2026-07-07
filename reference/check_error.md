# Error Catching

This function is used to check for errors in the
[`gformula`](https://adayim.github.io/causalMed/reference/gformula.md).

## Usage

``` r
check_error(data, id_var, base_vars, exposure, time_var, models)
```

## Arguments

- data:

  A `data.frame` (long format): one row per `id_var` per `time_var`.

- id_var:

  Character scalar. Subject identifier column name.

- base_vars:

  Character vector of time-fixed baseline covariates (may be empty).

- exposure:

  Character scalar. Exposure variable to intervene on (must be in
  `data`).

- time_var:

  Character scalar. Time variable column name (ordered;
  integer/numeric).

- models:

  A list of model specifications evaluated in temporal order. The order
  appeared in the list should reflect the temporal ordering of the
  variables, in another way data generation process. See
  [`spec_model`](https://adayim.github.io/causalMed/reference/spec_model.md)
  for a recommended constructor.

## Value

No value is returned.
