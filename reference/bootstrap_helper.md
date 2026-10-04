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
  future_seed = TRUE,
  time_seq
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

- mc_sample:

  Number of subjects in the simulated Monte Carlo cohort. The default,
  `NULL`, uses 50 times the number of subjects in `data` and reports the
  value unless `quiet = TRUE`. A larger value reduces Monte Carlo error,
  which falls as `1/sqrt(mc_sample)`; it does not change the estimand.

- R:

  Number of bootstrap replicates (default `500`); `R = 1` skips the
  bootstrap. Replicates run through
  [`future.apply::future_lapply`](https://future.apply.futureverse.org/reference/future_lapply.html):
  sequentially unless a parallel plan is set with
  [`future::plan()`](https://future.futureverse.org/reference/plan.html),
  e.g. `plan(multisession)`.

- future_seed:

  Logical or integer. Seed is passed to future_lapply.

- time_seq:

  Sorted distinct time points from the input data (required); see
  `.run_interventions`. Every replicate simulates this grid, so a
  resample missing a rare time value cannot shorten it.
