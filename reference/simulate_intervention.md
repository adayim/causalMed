# Monte Carlo simulation

Internal use only. Monte Carlo simulation.

## Usage

``` r
simulate_intervention(
  data,
  models,
  exposure,
  time_var,
  time_seq,
  intervention = NULL,
  init_recode = NULL,
  in_recode = NULL,
  out_recode = NULL,
  mediation_type = c(NA, "N", "I"),
  return_data = FALSE,
  med_pool = NULL,
  collect_pool = FALSE
)
```

## Arguments

- data:

  A `data.frame` in long format.

- models:

  List of
  [`spec_model`](https://adayim.github.io/causalMed/reference/spec_model.md)
  objects, in the order the variables are generated.

- exposure:

  Character. Name of the exposure to intervene on.

- time_var:

  Character. Name of the numeric time variable. Each distinct value is
  one simulated step, in ascending order.

- time_seq:

  Time sequence vector of the data.

- intervention:

  Named list of interventions. Each element is `NULL` (the natural
  course: exposure drawn from its fitted model), a 0/1 value or vector
  with one value per distinct time point (a static intervention), or a
  [`dyn_int`](https://adayim.github.io/causalMed/reference/dyn_int.md)
  rule, e.g.
  `list(natural = NULL, treat_if_high = dyn_int(as.numeric(L1 > 0)))`.
  `NULL` (the default) runs the natural course only.

- init_recode:

  [`recodes`](https://adayim.github.io/causalMed/reference/recodes.md)
  applied at the first time step, before the models are evaluated, e.g.
  to set lags to their baseline value.

- in_recode:

  [`recodes`](https://adayim.github.io/causalMed/reference/recodes.md)
  applied at the start of each later time step, before the models are
  evaluated, e.g. to update lags.

- out_recode:

  [`recodes`](https://adayim.github.io/causalMed/reference/recodes.md)
  applied at the end of each time step, the first included, e.g. to
  advance cumulative counts or carry an absorbing state forward.

- return_data:

  Logical. Return the simulated data (default `FALSE`; can be large).

- med_pool:

  Optional named list keyed by mediator response variable. Each element
  is itself a pre-permuted list of `T` vectors (one per time point, each
  of length `nrow(data)`, in the mediator's own type) supplying the
  cross-regime joint trajectory for that mediator. The per-time slice is
  delivered to `simulate_data`.

- collect_pool:

  Logical. If `TRUE`, capture the simulated trajectory of every mediator
  into a named list of per-time vector lists and return it alongside the
  risk estimate.
