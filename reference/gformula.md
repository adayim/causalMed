# Parametric g-formula for time-varying interventions

Estimates the mean outcome, or the risk for a survival outcome, that
would be seen if everyone followed each of several exposure strategies,
and the contrasts between them (the total effect), using the parametric
g-formula (Robins 1986). The models in `models` are fitted to the
observed data and then used to simulate each subject forward in time
under every intervention.

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

- intervention:

  Named list of interventions. Each element is `NULL` (the natural
  course: exposure drawn from its fitted model), a 0/1 value or vector
  with one value per distinct time point (a static intervention), or a
  [`dyn_int`](https://adayim.github.io/causalMed/reference/dyn_int.md)
  rule, e.g.
  `list(natural = NULL, treat_if_high = dyn_int(as.numeric(L1 > 0)))`.
  `NULL` (the default) runs the natural course only.

- ref_int:

  Reference for the contrasts: `0` or `"natural"` (default) for the
  natural course, or the position or name of an element of
  `intervention`. See Details.

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

- return_fitted:

  Logical. Return the full fitted model objects (default `FALSE`: calls
  and coefficients only).

- mc_sample:

  Number of subjects in the simulated Monte Carlo cohort. The default,
  `NULL`, uses 50 times the number of subjects in `data` and reports the
  value unless `quiet = TRUE`. A larger value reduces Monte Carlo error,
  which falls as `1/sqrt(mc_sample)`; it does not change the estimand.

- return_data:

  Logical. Return the simulated data (default `FALSE`; can be large).

- R:

  Number of bootstrap replicates (default `500`); `R = 1` skips the
  bootstrap. Replicates run through
  [`future.apply::future_lapply`](https://future.apply.futureverse.org/reference/future_lapply.html):
  sequentially unless a parallel plan is set with
  [`future::plan()`](https://future.futureverse.org/reference/plan.html),
  e.g. `plan(multisession)`.

- quiet:

  Logical. Suppress progress messages (default `FALSE`).

- seed:

  Integer random seed (default `12345`) for the Monte Carlo simulation
  and the bootstrap replicates; the global RNG state is restored on
  exit. `NULL` disables seeding, so repeated calls differ.

## Value

An object of class `"gformula"`, printed by
[`print.gformula`](https://adayim.github.io/causalMed/reference/print.gformula.md),
with components:

- `effect_size`: mean outcome (risk, for survival) under each
  intervention (`Intervention`, `Est`).

- `estimate`: with two or more interventions, the risk difference and
  risk ratio of each against `ref_int` (`Intervention`, `Risk_type`,
  `Estimate`).

- With `R > 1`, both tables gain `Sd` and percentile (`perct_lcl`,
  `perct_ucl`) and normal-approximation (`norm_lcl`, `norm_ucl`) limits,
  and `boot_estimates` holds the per-replicate estimates
  (`$interventions`, `$contrasts`).

- `sim_data`: with `return_data = TRUE`, the simulated data as an
  end-of-follow-up snapshot: one row per Monte Carlo subject per
  intervention, each variable at its last simulated time step, with the
  accumulated `Pred_Y`. Earlier time steps are not kept.

- `fitted_models`: the fitted models, as full objects when
  `return_fitted = TRUE`, otherwise their calls and coefficients.

- `data_summary`: numbers of subjects, rows and time points in `data`.

- `observed`: a nonparametric benchmark printed beside the estimates:
  the observed mean outcome at the last time point, or the product-limit
  cumulative incidence for a survival outcome.

- `call`, `all.args`: the matched call and the evaluated arguments.

## Details

**Data.** Long format: one row per subject per time point. For a
survival outcome the data must be in risk-set form, with every row after
the event and after loss to follow-up removed. A `"censor"` model may be
declared; it is simulated in the natural course but does not change any
reported risk, each of which is the risk under eliminated loss to
follow-up.

**Models.** `models` is a list of
[`spec_model`](https://adayim.github.io/causalMed/reference/spec_model.md)
objects in the order the variables are generated within a time point.
The `"outcome"` or `"survival"` model gives the predicted outcome
(`Pred_Y`) under each intervention: the fitted mean or risk of an
end-of-follow-up outcome, or the cumulative risk built from the fitted
hazards.

**Interventions.** `intervention` is a named list, for example
`list(natural = NULL, always = 1, never = 0)`. `ref_int = 0` or
`"natural"` (the default) compares against the natural course: the
`NULL` element if there is one, otherwise a `natural` element is added.
Because the natural course draws the exposure from its fitted model,
`models` then needs an `"exposure"` model. An integer position or a name
selects one of your interventions instead, and no natural course is
added.

**Rank-deficient models.** Terms that cannot be estimated (a collinear
term, or one that is constant in the rows used for fitting, such as
`time` in an outcome recorded only at the last time point) are dropped
from the simulation, as
[`predict`](https://rdrr.io/r/stats/predict.html) would drop them, and
named in the warning summary printed on exit.

The results rest on the correct temporal ordering and specification of
the models and on the usual g-formula assumptions (consistency,
positivity, no unmeasured confounding), which the package cannot check.

## References

Robins, J. M. (1986). A new approach to causal inference in mortality
studies with a sustained exposure period—application to control of the
healthy worker survivor effect. *Mathematical Modelling*, 7(9–12),
1393–1512.
[doi:10.1016/0270-0255(86)90088-6](https://doi.org/10.1016/0270-0255%2886%2990088-6)

Keil, A. P., Edwards, J. K., Richardson, D. B., Naimi, A. I., & Cole, S.
R. (2014). The parametric g-formula for time-to-event data: intuition
and a worked example. *Epidemiology*, 25(6), 889–897.
[doi:10.1097/EDE.0000000000000160](https://doi.org/10.1097/EDE.0000000000000160)

## See also

[`spec_model`](https://adayim.github.io/causalMed/reference/spec_model.md),
[`recodes`](https://adayim.github.io/causalMed/reference/recodes.md),
[`dyn_int`](https://adayim.github.io/causalMed/reference/dyn_int.md),
[`mediation`](https://adayim.github.io/causalMed/reference/mediation.md)
for direct and indirect effects, and
[`vignette("causalMed-03-gformula")`](https://adayim.github.io/causalMed/articles/causalMed-03-gformula.md).

## Examples

``` r
data(nonsurvivaldata)

# Models in the order the variables are generated: A -> L1 -> L2 -> Y
models <- list(
  spec_model(A ~ V + lag1_A + lag1_L1 + lag1_L2 + time,
             var_type = "binary", mod_type = "exposure"),
  spec_model(L1 ~ V + A + lag1_L1 + time,
             var_type = "normal", mod_type = "covariate"),
  spec_model(L2 ~ V + A + lag1_L2 + time,
             var_type = "binary", mod_type = "covariate"),
  spec_model(Y_bin ~ V + A + L1 + L2,
             var_type = "binary", mod_type = "outcome")
)

fit <- gformula(
  data = nonsurvivaldata, id_var = "id", time_var = "time",
  base_vars = "V", exposure = "A", models = models,
  intervention = list(natural = NULL, always = 1, never = 0),
  init_recode = recodes(lag1_A = 0, lag1_L1 = 0, lag1_L2 = 0),
  in_recode   = recodes(lag1_A = A, lag1_L1 = L1, lag1_L2 = L2),
  mc_sample = 2000,
  R = 1,          # use R > 1 (e.g. 500) for bootstrap confidence intervals
  quiet = TRUE
)
fit
#> Call:
#> gformula(data = nonsurvivaldata, id_var = "id", base_vars = "V", 
#>     exposure = "A", time_var = "time", models = models, intervention = list(natural = NULL, 
#>         always = 1, never = 0), init_recode = recodes(lag1_A = 0, 
#>         lag1_L1 = 0, lag1_L2 = 0), in_recode = recodes(lag1_A = A, 
#>         lag1_L1 = L1, lag1_L2 = L2), mc_sample = 2000, R = 1, 
#>     quiet = TRUE)
#> 
#> --- Analysis setup ---
#>   Exposure     : A
#>   Outcome      : Y_bin  [mean outcome at t = 4, end of follow-up]
#>   Time variable: time  (5 time points: 0 ... 4)
#>   ID variable  : id
#>   Baseline vars: V
#>   Data         : 3,000 individuals, 15,000 observations
#>   MC sample    : 2000
#>   Bootstrap R  : none
#>   Seed         : 12345
#>   Reference    : natural
#> 
#> --- Mean outcome by intervention --- 
#>    Intervention    Est
#>          <fctr>  <num>
#> 1:      natural 0.2258
#> 2:       always 0.2443
#> 3:        never 0.1022
#>   Observed (nonparametric) mean of Y_bin at t = 4 (end of follow-up): 0.2333
#>   (informal model check: compare with the natural-course intervention)
#> 
#> --- Contrasts vs. reference intervention --- 
#>        Intervention  Risk_type Estimate
#>              <char>     <char>    <num>
#> 1: always - natural Difference   0.0185
#> 2: always / natural      Ratio   1.0820
#> 3:  never - natural Difference  -0.1236
#> 4:  never / natural      Ratio   0.4525
```
