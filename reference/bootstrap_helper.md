# Calculate mediation analysis confidence interval

Used to calculate confidence interval using non-parametric bootstrap
methods.

## Usage

``` r
bootstrap_helper(
  data,
  id_var,
  base_vars,
  time_var,
  exposure,
  models,
  intervention,
  init_recode = NULL,
  in_recode = NULL,
  out_recode = NULL,
  mc_sample = 10000,
  mediation_type = c(NA, "N", "I"),
  n_vw = 1L,
  R = 500,
  progress_bar = TRUE,
  future_seed = TRUE
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

  If `intervention` is `NULL`, only the natural course is evaluated. A
  `natural` element is also added automatically when `ref_int` asks for
  the natural course and the list contains no `NULL` element; see
  `ref_int` and Details.

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

- mc_sample:

  Integer. Size of the Monte Carlo cohort simulated under each
  intervention, counted in subjects. Defaults to `NULL`, which resolves
  to 50 times the number of subjects in `data` and reports the value it
  chose unless `quiet = TRUE`. Monte Carlo error falls as
  `1/sqrt(mc_sample)` while sampling error falls with the number of
  subjects, so the two stay in a fixed ratio when `mc_sample` is set as
  a multiple of the subject count; a larger multiple buys precision in
  the point estimate and does not change what is being estimated.

- R:

  Number of bootstrap replicates. If `R > 1`, computation uses
  [`future.apply::future_lapply`](https://future.apply.futureverse.org/reference/future_lapply.html)
  and runs sequentially unless a parallel plan is set; see
  [`future::plan`](https://future.futureverse.org/reference/plan.html).
  Use `plan(multisession)` on Windows and `plan(multicore)` on
  Unix-alikes to enable parallel bootstrap. Default `500`.

- future_seed:

  Logical or integer. Seed is passed to future_lapply.
