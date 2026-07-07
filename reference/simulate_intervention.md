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

  A `data.frame` (long format): one row per `id_var` per `time_var`.

- models:

  A list of model specifications evaluated in temporal order. The order
  appeared in the list should reflect the temporal ordering of the
  variables, in another way data generation process. See
  [`spec_model`](https://adayim.github.io/causalMed/reference/spec_model.md)
  for a recommended constructor.

- exposure:

  Character scalar. Exposure variable to intervene on (must be in
  `data`).

- time_var:

  Character scalar. Time variable column name (ordered;
  integer/numeric).

- time_seq:

  Time sequence vector of the data.

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

- init_recode:

  Optional expression/function applied once at time 0 before the Monte
  Carlo loop (e.g., initializing baseline-derived variables). Should be
  defined with
  [`recodes`](https://adayim.github.io/causalMed/reference/recodes.md).
  See Details.

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

- return_data:

  Logical. If `TRUE`, return the stacked simulated data (all
  interventions) including predicted outcomes; may be large. Default
  `FALSE`.

- med_pool:

  Optional named list keyed by mediator response variable. Each element
  is a pre-permuted \`nrow(data) × T\` matrix (column names = time index
  as character) supplying the cross-world joint trajectory for that
  mediator. The per-time slice is delivered to `simulate_data`.

- collect_pool:

  Logical. If `TRUE`, capture the simulated trajectory of every mediator
  into a named list of \`nrow(data) × T\` matrices and return it
  alongside the risk estimate.
