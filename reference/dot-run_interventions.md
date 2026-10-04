# Main calculation function

This function will receive the parameters and fit model. After the model
is fitted, random samples will be drawn from the data and apply the
intervention.

## Usage

``` r
.run_interventions(
  data,
  id_var,
  base_vars,
  time_var,
  exposure,
  models,
  intervention,
  in_recode = NULL,
  out_recode = NULL,
  init_recode = NULL,
  mediation_type = c(NA, "N", "I"),
  mc_sample = 10000,
  n_vw = 1L,
  return_fitted = FALSE,
  return_data = FALSE,
  seed = NULL,
  time_seq = NULL
)
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

- time_var:

  Character. Name of the numeric time variable. Each distinct value is
  one simulated step, in ascending order.

- exposure:

  Character. Name of the exposure to intervene on.

- models:

  List of
  [`spec_model`](https://adayim.github.io/causalMed/reference/spec_model.md)
  objects, in the order the variables are generated.

- intervention:

  Named list of interventions. Each element is `NULL` (the natural
  course: exposure drawn from its fitted model), a 0/1 value or vector
  with one value per distinct time point (a static intervention), or a
  [`dyn_int`](https://adayim.github.io/causalMed/reference/dyn_int.md)
  rule, e.g.
  `list(natural = NULL, treat_if_high = dyn_int(as.numeric(L1 > 0)))`.
  `NULL` (the default) runs the natural course only.

- in_recode:

  [`recodes`](https://adayim.github.io/causalMed/reference/recodes.md)
  applied at the start of each later time step, before the models are
  evaluated, e.g. to update lags.

- out_recode:

  [`recodes`](https://adayim.github.io/causalMed/reference/recodes.md)
  applied at the end of each time step, the first included, e.g. to
  advance cumulative counts or carry an absorbing state forward.

- init_recode:

  [`recodes`](https://adayim.github.io/causalMed/reference/recodes.md)
  applied at the first time step, before the models are evaluated, e.g.
  to set lags to their baseline value.

- mediation_type:

  Type of the mediation analysis, if the value is `NA` no mediation
  analysis will be performed (default). It will be ignored if the
  intervention is not `NULL`

- mc_sample:

  Number of subjects in the simulated Monte Carlo cohort. The default,
  `NULL`, uses 50 times the number of subjects in `data` and reports the
  value unless `quiet = TRUE`. A larger value reduces Monte Carlo error,
  which falls as `1/sqrt(mc_sample)`; it does not change the estimand.

- n_vw:

  Integer. Number of independent permutation draws averaged for
  interventional pool-drawing interventions (Vansteelandt-Williamson
  repetition). Reference interventions (no mediator overrides) and
  natural-effect interventions are unaffected. Default `1L`;
  [`mediation()`](https://adayim.github.io/causalMed/reference/mediation.md)
  sets this to 2 to match the SAS mGFORMULA macro.

- return_fitted:

  Return the fitted model (default is FALSE).

- return_data:

  Logical. Return the simulated data (default `FALSE`; can be large).

- seed:

  Integer random seed (default `12345`) for the Monte Carlo simulation
  and the bootstrap replicates; the global RNG state is restored on
  exit. `NULL` disables seeding, so repeated calls differ.

- time_seq:

  Numeric vector: the sorted distinct time points to simulate, as
  computed by
  [`gformula()`](https://adayim.github.io/causalMed/reference/gformula.md)/[`mediation()`](https://adayim.github.io/causalMed/reference/mediation.md)
  from the *input* data. Passed to every bootstrap replicate so all
  passes simulate the same steps even when a resample lacks a rarely
  observed time value. When `NULL` the grid is derived from `data`.
