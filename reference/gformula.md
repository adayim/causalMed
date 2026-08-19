# Parametric G-formula for Time-varying Intervention Analyses

Implements the \*\*parametric g-formula\*\* to estimate counterfactual
mean outcomes under one or more user-specified exposure interventions in
longitudinal data. Supports settings with time-fixed baselines,
time-varying exposure/mediator/confounders, optional survival/terminal
outcomes, and nonparametric uncertainty via bootstrap.

## Usage

``` r
gformula(
  data,
  id_var,
  base_vars,
  exposure,
  time_var,
  models,
  intervention = NULL,
  ref_int = 0,
  init_recode = NULL,
  in_recode = NULL,
  out_recode = NULL,
  return_fitted = FALSE,
  mc_sample = NULL,
  return_data = FALSE,
  R = 500,
  quiet = FALSE,
  seed = 12345
)
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

- ref_int:

  Reference intervention for contrasts. Either an integer index (`0` =
  natural course; `1`, `2`, … = elements of `intervention`) or a
  character name matching an element (e.g., `"always"`). `0` and
  `"natural"` both resolve to the `NULL` element of `intervention` if
  there is one, and otherwise add a `natural` element — which requires
  an exposure model in `models`. See Details. Default: `0`.

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

- return_fitted:

  Logical. If `TRUE`, return full fitted model objects; otherwise, a
  light-weight summary (call and coefficients). Default `FALSE`.

- mc_sample:

  Integer. Size of the Monte Carlo cohort simulated under each
  intervention, counted in subjects. Defaults to `NULL`, which resolves
  to 50 times the number of subjects in `data` and reports the value it
  chose unless `quiet = TRUE`. Monte Carlo error falls as
  `1/sqrt(mc_sample)` while sampling error falls with the number of
  subjects, so the two stay in a fixed ratio when `mc_sample` is set as
  a multiple of the subject count; a larger multiple buys precision in
  the point estimate and does not change what is being estimated.

- return_data:

  Logical. If `TRUE`, return the stacked simulated data (all
  interventions) including predicted outcomes; may be large. Default
  `FALSE`.

- R:

  Number of bootstrap replicates. If `R > 1`, computation uses
  [`future.apply::future_lapply`](https://future.apply.futureverse.org/reference/future_lapply.html)
  and runs sequentially unless a parallel plan is set; see
  [`future::plan`](https://future.futureverse.org/reference/plan.html).
  Use `plan(multisession)` on Windows and `plan(multicore)` on
  Unix-alikes to enable parallel bootstrap. Default `500`.

- quiet:

  Logical. If `TRUE`, suppress progress messages/bars. Default `FALSE`.

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

## Value

An object of class `"gformula"` with components:

- `call`: the matched call.

- `all.args`: a named list of evaluated arguments for reproducibility.

- `effect_size`: `data.table` with columns `Intervention` and `Est`
  (intervention-specific mean outcome). If `R > 1`, also includes `Sd`,
  percentile CIs (`perct_lcl`, `perct_ucl`) and normal CIs (`norm_lcl`,
  `norm_ucl`).

- `estimate`: if multiple interventions are provided, a `data.table` of
  contrasts vs. `ref_int` (columns typically include `Intervention`,
  `Risk_type`, `Estimate`, and (if `R > 1`) CI columns).

- `sim_data`: if `return_data = TRUE`, the simulated Monte Carlo dataset
  stacked across interventions with an `Intervention` column (can be
  large). It is the **end-of-follow-up snapshot** — one row per Monte
  Carlo subject per intervention, holding each variable at its final
  simulated time step alongside the accumulated `Pred_Y` — not a
  row-per-time-point panel: the simulation overwrites each variable in
  place as it steps through time.

- `fitted_models`: a named list of fitted models. If
  `return_fitted = TRUE`, returns full model objects plus attributes
  (`recodes`, `subset`, `var_type`, `mod_type`); otherwise, a compact
  list with `call` and `coeff`.

- `boot_estimates`: when `R > 1`, a list of per-replicate bootstrap
  estimates: `$interventions` (columns `replicate`, `Intervention`,
  `Est`) and `$contrasts` (columns `replicate`, `Intervention`,
  `Risk_type`, `Estimate`; `NULL` with a single intervention). These are
  scalar summaries whose size is independent of the input data. `NULL`
  when `R <= 1`.

- `data_summary`: list with the number of individuals (`n_id`),
  observations (`n_obs`), and time points (`n_times`, `t_min`, `t_max`)
  of the input data.

- `observed`: list with the observed nonparametric benchmark of the
  outcome (`value`, `label`): the mean outcome at the last time point,
  or the product-limit cumulative incidence for survival outcomes.
  Printed alongside the simulated means as an informal model check.

## Details

The function evaluates a sequence of user-specified models (see
`models`) in \*\*temporal order\*\* to simulate counterfactual
trajectories via Monte Carlo, producing: (i) intervention-specific mean
outcomes, and (ii) contrasts vs. a reference intervention (risk
ratio/difference). When `R > 1`, percentile and normal-approximation
confidence intervals are computed from bootstrap resamples.

\*\*Data requirements\*\*

- Long format: one row per subject per time point.

- `id_var`: unique subject identifier; `time_var`: ordered time index.

- Final outcome must be well-defined at the last relevant time for each
  subject. For survival-like settings, rows after the event of interest
  should be removed (or a censoring model should be included and handled
  in `models`).

\*\*Interventions\*\* Provide a named list `intervention` with exposure
values per time (e.g.,
`list(natural = NULL, always = c(1,1,1), never = c(0,0,0))`). If
`intervention` is `NULL`, the function evaluates the natural course
only.

The reference is resolved before simulation, so that `ref_int` always
names an intervention that exists:

- `ref_int = 0` and `ref_int = "natural"` both ask for the natural
  course. The `NULL` element of `intervention` is used if one was
  supplied — whatever it is named — and otherwise a `natural` element is
  added to the list.

- An integer position, or a name matching an element, uses that
  intervention as the reference; no natural course is added.

Because the default is `ref_int = 0`, a list holding only static or
dynamic interventions still gains a natural course. The natural course
draws the exposure from its fitted model, so `models` must then include
a `mod_type = "exposure"` model; supplying one is required unless
`ref_int` names one of the interventions you provided.

\*\*Model specification\*\* Each element of `models` is typically
created by
[`spec_model`](https://adayim.github.io/causalMed/reference/spec_model.md)
and must include: (i) the model formula/call, (ii) a `mod_type`
indicating its role (`"exposure"`, `"covariate"`, `"outcome"`,
`"survival"`, or `"censor"`), and (iii) a `var_type` specifying the
variable type used for simulation/prediction (`"binary"`, `"normal"`,
`"categorical"`, and `"custom"`). The list order must reflect the
data-generating process (temporal ordering). The outcome model is
detected internally and used for computing predicted outcomes
(`Pred\_Y`) under each intervention.

\*\*Re-coding hooks\*\*

- `init_recode`: executed once at time 0 before simulation (initialize
  baselines).

- `in_recode`: executed at the start of each time step (e.g., entry-time
  logic).

- `out_recode`: executed at the end of each time step (e.g., cumulative
  counts, lags).

## Note

Final outcome should be consistently defined at the terminal time for
each subject. For survival-type applications, remove rows after the
event of interest (or include and model censoring appropriately). The
function may record warnings internally and print them on exit. Results
depend on correct temporal ordering, model specification, positivity,
and no unmeasured confounding assumptions customary for g-formula.

If a fitted model turns out to be rank deficient — a collinear term, or
a covariate that is constant among the rows actually used to fit it,
such as `time` in an outcome recorded only at the end of follow-up — the
unestimable terms are dropped from the simulation, exactly as
[`predict`](https://rdrr.io/r/stats/predict.html) would, and are named
in the warning summary printed on exit. Treat that warning as a prompt
to fix the formula.

## References

Robins, J. M. (1986). A new approach to causal inference in mortality
studies with a sustained exposure period—application to control of the
healthy worker survivor effect. *Mathematical Modelling*, 7(9–12),
1393–1512.

Keil, A. P., Edwards, J. K., Richardson, D. B., Naimi, A. I., & Cole, S.
R. (2014). The parametric g-formula for time-to-event data: intuition
and a worked example. *Epidemiology*, 25(6), 889–897.

## Examples

``` r
if (FALSE) { # \dontrun{
## Toy longitudinal data (long format)
data(nonsurvivaldata)

## Specify models in temporal order, e.g.:
mod1<-spec_model(A ~ A_lag1  +V,var_type= "binary",mod_type = "exposure")
mod2<-spec_model(L ~ A + L_lag1  +V,var_type  = "normal",mod_type = "covariate")
mod3<-spec_model(Y ~  A + L  + V,var_type = "binary",mod_type = "outcome")
models1<-list(mod1,mod2,mod3)

## Define interventions over T time points:
ints <- list(natural = NULL,
             always = 1
             )

fit <- gformula(
  data = sim_data0,
  id_var = 'id',
  base_vars = "V",
  exposure = "A",
  time_var = "time",
  models = models1,
  intervention = ints,
  ref_int = 1,
  init_recode = recodes(A_lag1 = 0, L_lag1 = 0),
  in_recode = recodes(A_lag1 = A, L_lag1 = L),
  out_recode = NULL,
  return_fitted = TRUE,
  mc_sample = 100000,
  return_data = TRUE,
  R = 500,
  quiet = TRUE,
  seed = 250817
)

print(fit)
summary(fit)
} # }
```
