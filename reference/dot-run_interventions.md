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
  seed = NULL
)
```

## Arguments

- data:

  A `data.frame` (long format): one row per `id_var` per `time_var`.

- id_var:

  Character scalar. Subject identifier column name.

- base_vars:

  Character vector of time-fixed baseline covariates (may be empty).

- time_var:

  Character scalar. Time variable column name (ordered;
  integer/numeric).

- exposure:

  Character scalar. Exposure variable to intervene on (must be in
  `data`).

- models:

  A list of model specifications evaluated in temporal order. The order
  appeared in the list should reflect the temporal ordering of the
  variables, in another way data generation process. See
  [`spec_model`](https://adayim.github.io/causalMed/reference/spec_model.md)
  for a recommended constructor.

- intervention:

  A named list specifying exposure interventions. Each element is one
  of:

  - `NULL` — the natural course (exposure drawn from its fitted model).

  - A numeric/logical scalar or vector (length 1 or equal to the number
    of time points) — a static intervention setting the exposure to that
    value at every (or each specific) time step.

  - A
    [`dyn_int`](https://adayim.github.io/causalMed/reference/dyn_int.md)
    object — a dynamic (rule-based) intervention whose expression is
    evaluated inside the simulated dataset at each time step. Column
    names (including the exposure after its natural-course draw) are in
    scope, e.g.
    `list(natural = NULL, threshold = dyn_int(as.numeric(A > 0)))`.

  If `intervention` is `NULL`, only the natural course is evaluated. If
  a `natural` intervention is not provided, it is added automatically
  and `ref_int` is set to `"natural"`.

- in_recode:

  Optional expression/function applied at the \*\*start\*\* of each time
  step (e.g., entry-time functional forms). Should be defined with
  [`recodes`](https://adayim.github.io/causalMed/reference/recodes.md).
  See Details.

- out_recode:

  Optional expression/function applied at the \*\*end\*\* of each time
  step (e.g., create lags, cumulative counts). Should be defined with
  [`recodes`](https://adayim.github.io/causalMed/reference/recodes.md).
  See Details.

- init_recode:

  Optional expression/function applied once at time 0 before the Monte
  Carlo loop (e.g., initializing baseline-derived variables). Should be
  defined with
  [`recodes`](https://adayim.github.io/causalMed/reference/recodes.md).
  See Details.

- mediation_type:

  Type of the mediation analysis, if the value is `NA` no mediation
  analysis will be performed (default). It will be ignored if the
  intervention is not `NULL`

- mc_sample:

  Integer. Monte Carlo sample size used for simulation. Default
  `nrow(data) * 100`.

- n_vw:

  Integer. Number of independent permutation draws averaged for
  interventional cross-world interventions (Vansteelandt-Williamson
  repetition). Reference interventions (no mediator overrides) and
  natural-effect interventions are unaffected. Default `1L`;
  [`mediation()`](https://adayim.github.io/causalMed/reference/mediation.md)
  sets this to 2 to match the SAS mGFORMULA macro.

- return_fitted:

  Return the fitted model (default is FALSE).

- return_data:

  Logical. If `TRUE`, return the stacked simulated data (all
  interventions) including predicted outcomes; may be large. Default
  `FALSE`.

- seed:

  Integer random seed for reproducibility. Default `12345`. The seed
  fixes the RNG stream used for the Monte Carlo simulation *and* the
  bootstrap replicates, and the caller's global RNG state is saved and
  restored on exit, so a seeded call does not disturb the ambient random
  stream. Pass `NULL` to disable seeding entirely: no
  [`set.seed()`](https://rdrr.io/r/base/Random.html) is called, the
  Monte Carlo draws consume and advance the ambient RNG stream, and
  repeated calls therefore give *different* results (reproducible only
  via an outer [`set.seed()`](https://rdrr.io/r/base/Random.html)). Use
  `seed = NULL` inside simulation loops that manage their own seeds.
