# G-formula based mediation analysis

Decompose the effect of a time-varying exposure into direct and indirect
components using the parametric mediational g-formula. Estimates
interventional effects (`"I"`, the default) or natural effects (`"N"`),
with bootstrap or influence-curve confidence intervals.

## Usage

``` r
mediation(
  data,
  id_var,
  base_vars,
  exposure,
  outcome,
  time_var,
  models,
  init_recode = NULL,
  in_recode = NULL,
  out_recode = NULL,
  mc_sample = NULL,
  mediation_type = c("I", "N"),
  n_vw = 2L,
  estimator = c("gcomp", "tmle"),
  tmle_weight_trunc = 0.995,
  return_fitted = FALSE,
  return_data = FALSE,
  R = 500,
  quiet = FALSE,
  seed = 12345L
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

- outcome:

  Character scalar. Name of the outcome variable in `data`. Must match
  the response variable in the outcome or survival model.

- time_var:

  Character scalar. Time variable column name (ordered;
  integer/numeric).

- models:

  A list of model specifications evaluated in temporal order. The order
  appeared in the list should reflect the temporal ordering of the
  variables, in another way data generation process. See
  [`spec_model`](https://adayim.github.io/causalMed/reference/spec_model.md)
  for a recommended constructor.

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

- mediation_type:

  Character. Type of mediation effect: `"I"` for interventional effect
  (default) or `"N"` for natural effect. **Identifiability.** Natural
  effects (`"N"`) are *not identifiable* from observational data when a
  time-varying confounder of the mediator-outcome relationship is itself
  affected by prior exposure (Avin, Shpitser & Pearl 2005; VanderWeele
  2014; VanderWeele & Tchetgen Tchetgen 2017). For that setting
  VanderWeele & Tchetgen Tchetgen (2017) propose the randomized
  interventional analogues, which `"I"` targets and which remain
  identifiable — hence the default. Choosing `"N"` when a covariate
  model carries the exposure on its right-hand side triggers a warning,
  repeated by [`print()`](https://rdrr.io/r/base/print.html). That check
  reads formulas only: it neither establishes that such a covariate
  confounds the mediator-outcome relationship nor, by its silence, that
  no intermediate confounder exists. Judging the causal structure
  remains the analyst's responsibility. **Interpretation of `"I"`.**
  Interventional indirect effects buy identifiability at a price: absent
  stronger assumptions they do not satisfy the sharp null criterion of
  Miles (2023) — being null whenever no individual-level indirect effect
  exists — which the natural indirect effect does satisfy. A non-zero
  IIE therefore does not by itself show that the mediator transmits the
  effect for any individual. The choice is substantive, not
  computational.

- n_vw:

  Integer. Number of independent permutation draws averaged for each
  intervention that draws its mediators from a permuted pool. Under
  `mediation_type = "I"` that is *every* intervention in the
  decomposition — the references `Phi00` and `Phi11` just as much as the
  cross-regime `Phi10` and `Phi1_k` — but not the fixed-exposure,
  natural-mediator interventions `nat0`/`nat1`, whose mediators come
  from their own fitted models and involve no permutation. The same
  averaging is applied within every bootstrap replicate. Default `2L` to
  match the SAS `mGFORMULA` macro's `n_vw = 2`. Set to `1L` to disable
  averaging (faster, slightly noisier Monte Carlo). Has no effect when
  `mediation_type = "N"` (the natural-effect mediator swap is not
  permutation-based).

- estimator:

  Character. `"gcomp"` (default) for the parametric g-formula plug-in
  (Monte Carlo simulation + bootstrap CIs), or `"tmle"` for the targeted
  minimum loss-based estimator of Zheng & van der Laan (2017, Section
  4.3). `"tmle"` is available for `mediation_type = "N"` only (single
  mediator, binary {0, 1} outcome or survival event indicator; an
  exposure model is required and `var_type = "custom"` /
  `spec_model(subset = )` are not supported). It uses backward iterated
  regressions with logistic fluctuations rather than forward simulation:
  multiply robust (consistent when specific subsets of the nuisance
  models are correct, not only when all are), Wald CIs from the
  efficient influence curve without bootstrapping (`R` ignored), and no
  Monte Carlo error (`mc_sample` ignored). **Recode restriction:** only
  lag-style recodes are evaluated. `in_recode` entries must copy a
  single column (e.g. `recodes(lag_A = A)`); exposure lags must copy the
  exposure itself (chained lags such as `recodes(lag2_A = lag_A)` are
  rejected); `init_recode` entries must be a constant or a column name;
  and `out_recode` is unsupported. Derived recodes (splines, cumulative
  counts, carry-forward flags) and deeper exposure history require
  `"gcomp"`. Violations error rather than pass silently. **Working-model
  form:** your formulas are used as written for the conditional
  densities behind the clever covariates, but the targeted sequential
  regressions are built as *additive main-effects* models in the
  variables those formulas name, so transformations and interactions
  ([`poly()`](https://rdrr.io/r/stats/poly.html), splines, `A:M`) are
  not carried into them. That is a property of this implementation, not
  of the published estimator, whose results assume correctly specified
  nuisance models. Positivity violations (few subjects following an
  intervened regime) skip the affected fluctuation steps and truncate
  extreme weights, with collected warnings — inspect those before
  trusting the affected functionals.

- tmle_weight_trunc:

  Numeric in (0, 1\]. Quantile at which each clever-covariate weight
  vector is truncated in the TMLE fluctuation steps (positivity
  protection). Default `0.995`; `1` disables truncation. Ignored for
  `estimator = "gcomp"`.

- return_fitted:

  Logical. If `TRUE`, return full fitted model objects; otherwise, a
  light-weight summary (call and coefficients). Default `FALSE`.

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

  Integer random seed. Default `12345L`. Fixes the stream for both the
  Monte Carlo simulation and the bootstrap replicates (the latter via
  future.apply's `future_seed`); the caller's global RNG state is
  restored on exit. Pass `NULL` to disable seeding — draws then consume
  the ambient stream and repeated calls differ — as needed inside
  simulation loops that manage their own seeds. Give concurrent analyses
  distinct seeds to keep their bootstrap draws independent.

## Value

An object of class `"gformula"` with components:

- `call`: the matched function call.

- `all.args`: a named list of evaluated arguments for reproducibility.

- `effect_size`: a `data.table` with one row per simulated intervention
  (columns `Intervention` and `Est`). For `mediation_type = "I"` the
  mediator(s) in every decomposition intervention are drawn from
  independently-permuted marginal pools (stochastic draws \\G\\):
  `Phi00` (= \\E\[Y\_{0,G_0}\]\\) and `Phi11` (= \\E\[Y\_{1,G_1}\]\\)
  are the interventional reference interventions, `Phi10` (=
  \\E\[Y\_{1,G_0}\]\\) is the cross-regime intervention, and for \\N \ge
  2\\ mediators `Phi1_k` (k = 1, …, N-1) capture the sequential
  per-mediator transition. Two additional fixed-exposure,
  natural-mediator interventions `nat0` (= \\E\[Y_0\]\\) and `nat1` (=
  \\E\[Y_1\]\\) give the plug-in total effect. For
  `mediation_type = "N"` the interventions `Phi00`/`Phi11` are the
  natural never-/always-treat interventions (no permutation). With
  `R > 1` the table also carries `Sd`, `perct_lcl`/`perct_ucl`, and
  `norm_lcl`/`norm_ucl`.

- `estimate`: a `data.table` summarizing the decomposition of effects.
  For a single mediator under `mediation_type = "I"` the rows are:

  - `"Indirect effect"` = \\E\[Y\_{1,G_1}\] - E\[Y\_{1,G_0}\]\\

  - `"Direct effect"` = \\E\[Y\_{1,G_0}\] - E\[Y\_{0,G_0}\]\\

  - `"Total effect"` = \\E\[Y_1\] - E\[Y_0\]\\ (natural plug-in
    g-formula; *not* the sum of the components)

  - `"TE - (Direct + Indirect)"` = the decomposition residual \\TE -
    (IDE + IIE)\\ = \\TE\\ minus the interventional overall effect
    \\E\[Y\_{1,G_1}\] - E\[Y\_{0,G_0}\]\\ (generally non-zero for
    interventional effects; this row is *absent* for
    `mediation_type = "N"`, where the decomposition sums exactly to
    \\TE\\).

  - `"Mediation Proportion"` = \\\sum_k IIE(M_k) / OE \times 100\\\\,
    where \\OE = IDE + \sum_k IIE(M_k)\\ is the interventional overall
    effect. Numerator and denominator come from the same decomposition,
    so \\IDE/OE\\ and this proportion sum to 100%. It is a share of the
    *overall* effect, **not** of the natural plug-in total effect: the
    two differ by the decomposition residual, which is reported in its
    own row on the RD scale. This is the quantity Lin et al. (2017,
    *Stat Med*, Table 2) report; their design has no separate
    fixed-exposure, natural-mediator intervention, so the "total effect"
    in their formula is \\\Phi\_{11}-\Phi\_{00}\\. Zheng & van der
    Laan (2017) do not define a proportion mediated; under
    `mediation_type = "N"` there is no residual and the two denominators
    coincide.

  A separate multiplicative proportion is *not* reported. The risk-ratio
  formula \\RR\_{IDE}(\prod_k RR\_{IIE_k}-1)/(RR\_{OE}-1)\\ is not a
  second scale: it simplifies to
  \\(\Phi\_{11}-\Phi\_{10})/(\Phi\_{11}-\Phi\_{00})\\, the same number
  as the proportion above. A genuine ratio-scale proportion would need a
  scale-specific definition.

  For multiple mediators each indirect effect is labelled
  `"Indirect effect (<mediator>)"` and the proportion sums them
  (Yamamuro et al. 2021). Columns: `RD`, `RR`; with `R > 1` also `Sd`,
  percentile CIs (`perct_lcl`/`perct_ucl`), and normal CIs
  (`norm_lcl`/`norm_ucl`) for RD; and `Sd_RR`,
  `perct_lcl_RR`/`perct_ucl_RR`, `norm_lcl_RR`/`norm_ucl_RR` for RR.
  `RR` is `NA` for the proportion rows.

- `sim_data`: if `return_data = TRUE`, the simulated Monte Carlo dataset
  used internally (can be large), stacked across interventions with an
  `Intervention` column. It is the **end-of-follow-up snapshot** — one
  row per Monte Carlo subject per intervention, holding each variable at
  its final simulated time step alongside the accumulated `Pred_Y` — not
  a row-per-time-point panel. With `n_vw > 1` the pool-drawing
  interventions are simulated `n_vw` times and only the *last*
  permutation is retained here, whereas `effect_size$Est` averages all
  `n_vw` of them. Recomputing `mean(Pred_Y)` from `sim_data` will
  therefore not reproduce `Est` exactly for those interventions (it does
  for `nat0`/`nat1`); use `n_vw = 1` if you need the two to agree.

- `fitted_models`: a named list of fitted models keyed by outcome,
  exposure, and mediator variables. If `return_fitted = TRUE`, returns
  full model objects plus attributes (`recodes`, `subset`, `var_type`,
  `mod_type`); otherwise, a compact list with `call` and `coeff`.

- `boot_estimates`: when `R > 1`, a list of per-replicate bootstrap
  estimates: `$interventions` (columns `replicate`, `Intervention`,
  `Est`) and `$effects` (columns `replicate`, `Effect`, `RD`, `RR`;
  includes the per-replicate Mediation Proportion draws). Useful for
  diagnostics such as counting non-finite replicates or computing custom
  intervals. These are scalar summaries whose size is independent of the
  input data. `NULL` when `R <= 1`.

- `data_summary`: list with the number of individuals (`n_id`),
  observations (`n_obs`), and time points (`n_times`, `t_min`, `t_max`)
  of the input data.

- `observed`: list with the observed nonparametric benchmark of the
  outcome (`value`, `label`): the mean outcome at the last time point,
  or the product-limit cumulative incidence for survival outcomes.
  Printed next to the simulated intervention means as an informal
  benchmark.

- `tmle_diag`: for `estimator = "tmle"` only, a list with the
  subject-level efficient-influence-curve matrix (`eic`, one column per
  functional), its column means (`eic_mean`, near zero for
  well-supported functionals at the targeted fit), and the number of
  subjects (`n`). `NULL` for `estimator = "gcomp"`.

- `intermediate_confounders`: for `mediation_type = "N"`, a character
  vector of covariate names whose model includes the exposure
  (exposure-affected/intermediate confounders that make the natural
  effects non-identifiable); empty when none. The `print` method
  re-surfaces a short identifiability caveat when this is non-empty.

## Details

\*\*Data requirements\*\* The input dataset must be in longitudinal
("long") format: one record per subject per time period. In survival
settings, all records after the event must be removed.

The exposure variable must be \*\*binary, coded as 0
(reference/untreated) and 1 (active/treated)\*\*. Methods implemented
here (Lin et al. 2017; Zheng & van der Laan 2017) are defined for binary
exposures only.

\*\*Model specification\*\* Each element of `models` is created by
[`spec_model`](https://adayim.github.io/causalMed/reference/spec_model.md)
and must specify: (i) a formula, (ii) `mod_type` (`"exposure"`,
`"covariate"`, `"mediator"`, `"outcome"`, `"survival"`, or `"censor"`),
and (iii) `var_type` (`"binary"`, `"normal"`, `"categorical"`, or
`"custom"`). At least one `"mediator"` model is required. Multiple
mediator models are supported under `mediation_type = "I"` (Yamamuro et
al. 2021); list them in temporal order. Multi-mediator analyses are not
supported under `mediation_type = "N"`.

A `"censor"` model declares a discrete-time censoring process. Under
every intervention the censoring indicator is set to zero, so the
estimand is the risk under eliminated loss to follow-up. With the
default `estimator = "gcomp"` the censoring model does not alter any
reported risk: no simulated subject is removed on any pass, and the risk
is the average of the analytic \\1-\prod_t(1-\hat p_t)\\ rather than of
the simulated event indicator. Its presence marks the analysis as a
survival setting, and it is used by the targeted estimator
(`estimator = "tmle"`).

The list order determines the simulation sequence at each time step and
must match your assumed data-generating process. A common ordering is
**A(t) -\> L(t) -\> M(t) -\> S(t)** (confounders not affected by
mediator) or **A(t) -\> M(t) -\> L(t) -\> S(t)** (confounders affected
by both exposure and mediator). The function checks that the mediator
precedes the outcome and that exposure precedes the mediator, and warns
if these are violated.

For non-standard covariate distributions (bounded, zero-inflated,
truncated), use `var_type = "custom"` with `custom_fit` and `custom_sim`
arguments to
[`spec_model`](https://adayim.github.io/causalMed/reference/spec_model.md).
See the custom covariate distributions section of
[`vignette("causalMed-03-gformula")`](https://adayim.github.io/causalMed/articles/causalMed-03-gformula.md).

\*\*Mediator pool (interventional effects)\*\* For
`mediation_type = "I"`, one pass at each treatment level \\a^\*\\ —
exposure held fixed, mediator following its fitted model — stores every
simulated individual's full trajectory \\M(1{:}T)\\ in a pool. Each
decomposition intervention that fixes a mediator to its \\a^\*\\ value,
*including the references* `Phi00` (\\a^\* = 0\\) and `Phi11` (\\a^\* =
1\\), permutes that pool once and assigns subject \\i\\ the whole
trajectory of pool individual \\\pi(i)\\ — a joint stochastic draw
\\G\_{a^\*}\\. Mediators are permuted independently. Every intervention
in the decomposition therefore uses the randomized mediator
distribution, not only the cross-regime ones.

The permutation is uniform over the whole simulated cohort, so the draw
is marginal over the baseline covariates as well as over \\i\\'s own
confounder history. VanderWeele and Tchetgen Tchetgen (2017) define both
a within-stratum draw \\G\_{a\|v}\\ and a whole-population draw
\\G\_{a}\\ with its own identifying formula; this package targets the
second, as does the fixed-endpoint algorithm of Lin et al. (2017,
*Epidemiology*), of which averaging over several permutations (`n_vw`)
is a step. Yamamuro et al. (2021, Eq. 2) and Lin et al. (2017, *Stat
Med*, Eq. 4) are written for the conditional form, but the algorithms in
both papers — and the SAS macros distributed with them — permute the
whole cohort. The forms agree when the counterfactual mediator
distribution does not depend on the baseline covariates, or when the
outcome mean is additively separable in them and the mediator; otherwise
they are different functionals.

**Survival outcomes.** This package implements the estimation algorithm
of Lin et al. (2017, *Stat Med*, Section 4): the mediator models are
fitted among those at risk, no one is then removed from the simulated
cohort, cumulative risk is accumulated analytically over the
discrete-time hazards, and the pool is permuted with equal weight. The
natural plug-in total effect is obtained from the separate
fixed-exposure, natural-mediator interventions (`nat0`, `nat1`), so
\\TE\\ need not equal the sum of the interventional direct and indirect
effects; their difference is reported as the decomposition residual.

\*\*Warnings\*\* Warnings from model fitting (e.g., convergence,
near-separation) are collected during the run and printed as a
deduplicated summary at function exit, including a repeat count when the
same message fires across bootstrap replicates.

\*\*Re-coding hooks\*\*

- `init_recode`: executed once at time 0 before simulation (initialize
  baselines).

- `in_recode`: executed at the start of each time step (e.g., entry-time
  logic).

- `out_recode`: executed at the end of each time step (e.g., cumulative
  counts, lags).

## References

Lin, S. H., Young, J. G., Logan, R., & VanderWeele, T. J. (2017).
Mediation analysis for a survival outcome with time-varying exposures,
mediators, and confounders. *Statistics in Medicine*, 36(26), 4153–4166.
[doi:10.1002/sim.7426](https://doi.org/10.1002/sim.7426)

VanderWeele, T. J., & Tchetgen Tchetgen, E. J. (2017). Mediation
analysis with time varying exposures and mediators. *Journal of the
Royal Statistical Society: Series B*, 79(3), 917–938.
[doi:10.1111/rssb.12194](https://doi.org/10.1111/rssb.12194)

Yamamuro, S., Shinozaki, T., Iimuro, S., & Matsuyama, Y. (2021).
Mediational g-formula for time-varying treatment and repeated-measured
multiple mediators: Application to atorvastatin's effect on
cardiovascular disease via cholesterol lowering and anti-inflammatory
actions in elderly type 2 diabetics. *Statistical Methods in Medical
Research*, 30(8), 1782–1799.
[doi:10.1177/09622802211025988](https://doi.org/10.1177/09622802211025988)

Zheng, W., & van der Laan, M. (2017). Longitudinal mediation analysis
with time-varying mediators and exposures, with application to survival
outcomes. *Journal of Causal Inference*, 5(2).
[doi:10.1515/jci-2016-0006](https://doi.org/10.1515/jci-2016-0006)

Miles, C. H. (2023). On the causal interpretation of randomised
interventional indirect effects. *Journal of the Royal Statistical
Society: Series B*, 85(4), 1154–1172.
[doi:10.1093/jrsssb/qkad066](https://doi.org/10.1093/jrsssb/qkad066)

## Examples

``` r
if (FALSE) { # \dontrun{
data(nonsurvivaldata)

models <- list(
  cov_model = spec_model(L ~ V+L_lag1+A+time,var_type= "normal",mod_type = "covariate"),
  mediator_model = spec_model(M ~ V + A + L + M_lag1 + time,
                              var_type = "normal", mod_type = "mediator"),
  outcome_model = spec_model(Y ~ V+A+M+A * M+L,  var_type= "binary",mod_type ="outcome")
  )

fit <- mediation(
 data = nonsurvivaldata,
 id_var = "id",
 base_vars = "V",
 exposure = "A",
 outcome = "Y",
 time_var = "time",
 models = models,
 init_recode = recodes(M_lag1=0,L_lag1=0),
 in_recode = recodes(M_lag1=M,L_lag1=L),
 mediation_type = "I",
 mc_sample = 100000,
 R = 500,
 return_data = FALSE,
 return_fitted = FALSE
 )

print(fit)
} # }
```
