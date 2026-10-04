# Mediation analysis with the parametric g-formula

Decomposes the effect of a time-varying binary exposure into a direct
effect and an indirect effect through one or more time-varying
mediators, using the parametric mediational g-formula. Two estimands are
available: interventional effects (`mediation_type = "I"`, the default;
Lin et al. 2017, VanderWeele & Tchetgen Tchetgen 2017) and natural
effects (`"N"`; Zheng & van der Laan 2017).

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
  exposure_regime = 1,
  reference_regime = 0,
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

- outcome:

  Character. Name of the outcome variable; must be the response of the
  `"outcome"` or `"survival"` model.

- time_var:

  Character. Name of the numeric time variable. Each distinct value is
  one simulated step, in ascending order.

- models:

  List of
  [`spec_model`](https://adayim.github.io/causalMed/reference/spec_model.md)
  objects, in the order the variables are generated.

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

- mediation_type:

  `"I"` (default) for interventional effects or `"N"` for natural
  effects.

  Zheng & van der Laan (2017, Lemma 1) identify the `"N"` effects, which
  draw the mediator from its conditional distribution under the other
  regime given each subject's own history, under sequential
  randomization and positivity. Reading them as *individual-level*
  natural effects additionally requires a cross-world independence
  assumption that is not expected to hold when a confounder of the
  mediator-outcome relationship is affected by prior exposure (Avin,
  Shpitser & Pearl 2005; VanderWeele & Tchetgen Tchetgen 2017), the
  setting for which VanderWeele & Tchetgen Tchetgen (2017) propose the
  interventional effects of `"I"`. With `"N"`, a warning (repeated by
  [`print()`](https://rdrr.io/r/base/print.html)) names any covariate
  whose model includes the exposure. The check reads formulas only: it
  does not establish that such a covariate confounds the
  mediator-outcome relationship, and its silence does not establish that
  none does.

  Interventional indirect effects do not, absent stronger assumptions,
  satisfy the sharp null criterion of Miles (2023), which the natural
  indirect effect does: a non-zero interventional indirect effect does
  not by itself show that the mediator transmits the effect for any
  individual.

- exposure_regime:

  Exposure regime \\a(1{:}T)\\ whose effect is decomposed: 0/1 values,
  either one value for every time point or one per distinct value of
  `time_var`, in sorted order. Default `1` (always exposed).

- reference_regime:

  Reference regime \\a^\*(1{:}T)\\, in the same form. Default `0` (never
  exposed). Must differ from `exposure_regime`.

- n_vw:

  Integer. Number of mediator permutations averaged for each
  pool-drawing intervention under `"I"` (`Phi00`, `Phi10`, `Phi11`,
  `Phi1_k`; not `nat0`/`nat1`), within every bootstrap replicate too.
  Default `2L`, as in the SAS `mGFORMULA` macro; `1L` is faster with
  more Monte Carlo noise. Ignored under `"N"`.

- estimator:

  `"gcomp"` (default): the parametric g-formula plug-in with bootstrap
  intervals. `"tmle"`: the targeted minimum loss-based estimator of
  Zheng & van der Laan (2017, Section 4.3), with Wald intervals from the
  efficient influence curve (`R` and `mc_sample` are ignored).

  `"tmle"` requires `mediation_type = "N"`, the default regimes, a
  single mediator, a binary outcome or survival indicator, and an
  exposure model. It accepts only lag-style recodes: `in_recode` entries
  that copy one column (an exposure lag must copy the exposure itself),
  `init_recode` entries that are a constant or a column name, no
  `out_recode`, no `subset`, no `var_type = "custom"`, and no model
  `recode` that reads the exposure or its lags. Other inputs raise an
  error. Its targeted regressions are additive main-effects working
  models in the variables your formulas name (transformations and
  interactions such as [`poly()`](https://rdrr.io/r/stats/poly.html) or
  `A:M` are not carried over), so the multiple robustness Zheng & van
  der Laan establish for correctly specified nuisance models is not
  claimed for this implementation. Where few subjects follow a regime,
  the affected fluctuation steps are skipped and weights truncated, with
  a warning.

- tmle_weight_trunc:

  Quantile in (0, 1\] at which the TMLE clever-covariate weights are
  truncated. Default `0.995`; `1` disables truncation. Ignored unless
  `estimator = "tmle"`.

- return_fitted:

  Logical. Return the full fitted model objects (default `FALSE`: calls
  and coefficients only).

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

  Integer random seed (default `12345L`) for the Monte Carlo simulation
  and the bootstrap replicates; the global RNG state is restored on
  exit. `NULL` disables seeding, so repeated calls differ. Give
  concurrent analyses different seeds.

## Value

An object of class `"gformula"`, printed by
[`print.gformula`](https://adayim.github.io/causalMed/reference/print.gformula.md),
with components:

- `effect_size`: mean outcome (risk, for survival) under each
  intervention (`Intervention`, `Est`). Under `"I"`: `Phi00` \\=
  E\[Y\_{0,G_0}\]\\, `Phi10` \\= E\[Y\_{1,G_0}\]\\, `Phi11` \\=
  E\[Y\_{1,G_1}\]\\ (mediator from a permuted pool), and `Phi1_k` (k =
  1, ..., N-1) with N mediators; under `"N"`, `Phi00` and `Phi11` use
  each subject's own mediator. Always `nat0` \\= E\[Y_0\]\\ and `nat1`
  \\= E\[Y_1\]\\. With non-default regimes, 0 and 1 in these labels
  stand for \\a^\*\\ and \\a\\.

- `estimate`: the decomposition, on the risk-difference (`RD`) and
  risk-ratio (`RR`) scales. Rows: `"Indirect effect"` (`Phi11 - Phi10`;
  one row per mediator, `"Indirect effect (<mediator>)"`, with several),
  `"Direct effect"` (`Phi10 - Phi00`), `"Total effect"` (`nat1 - nat0`),
  `"TE - (Direct + Indirect)"` (`"I"` only), and
  `"Mediation Proportion"`: the summed indirect effects as a percentage
  of the interventional overall effect `Phi11 - Phi00`, the quantity Lin
  et al. (2017, *Stat Med*, Table 2) report. It is not a share of the
  total effect, from which it differs by the residual. No separate
  risk-ratio proportion is given: \\RR\_{IDE}(\prod_k
  RR\_{IIE_k}-1)/(RR\_{OE}-1)\\ simplifies to the same number. `RR` is
  `NA` for the residual and proportion rows.

- With `R > 1`, both tables gain `Sd` and percentile (`perct_lcl`,
  `perct_ucl`) and normal-approximation (`norm_lcl`, `norm_ucl`) limits,
  with an `_RR` suffix on the risk-ratio scale; `boot_estimates` holds
  the per-replicate estimates (`$interventions`, `$effects`).

- `sim_data`: with `return_data = TRUE`, the simulated data as an
  end-of-follow-up snapshot: one row per Monte Carlo subject per
  intervention, each variable at its last simulated time step, with the
  accumulated `Pred_Y`. With `n_vw > 1` only the last permutation is
  kept while `effect_size$Est` averages all of them, so `mean(Pred_Y)`
  reproduces `Est` for the pool-drawing interventions only when
  `n_vw = 1` (always for `nat0`/`nat1`).

- `fitted_models`: the fitted models, as full objects when
  `return_fitted = TRUE`, otherwise their calls and coefficients.

- `data_summary`: numbers of subjects, rows and time points in `data`,
  and `regime_support`, the number of observed subjects whose exposure
  matches each regime at every time they were observed (`n_following`)
  and among those observed at every time point (`n_complete`). A subject
  with a missing exposure counts for neither regime. The estimation does
  not use these counts.

- `observed`: a nonparametric benchmark printed beside the estimates:
  the observed mean outcome at the last time point, or the product-limit
  cumulative incidence for a survival outcome.

- `intermediate_confounders`: under `"N"`, covariates whose model
  includes the exposure (see `mediation_type`).

- `tmle_diag`: with `estimator = "tmle"`, the subject-level efficient
  influence curve (`eic`), its column means (`eic_mean`) and the number
  of subjects (`n`).

- `call`, `all.args`: the matched call and the evaluated arguments.

## Details

**Data.** Long format: one row per subject per time point. For a
survival outcome, remove every row after the event and after loss to
follow-up. The exposure must be coded 0/1; the methods implemented here
are defined for binary exposures only.

**Models.** `models` is a list of
[`spec_model`](https://adayim.github.io/causalMed/reference/spec_model.md)
objects in the order the variables are generated within a time point,
for example exposure, confounders, mediator, outcome (**A(t) -\> L(t)
-\> M(t) -\> Y(t)**), or **A -\> M -\> L -\> Y** when the confounders
respond to the mediator. At least one `"mediator"` model is required.
Several mediators, listed in temporal order, are supported under `"I"`
only (Yamamuro et al. 2021). A warning is given if the exposure is
listed after the mediator or the mediator after the outcome.

**Exposure regimes.** The effects contrast two *static* regimes, \\a\\
(`exposure_regime`) and \\a^\*\\ (`reference_regime`). The defaults,
always versus never exposed, are the example Lin et al. (2017, *Stat
Med*, Section 2.2) and Zheng & van der Laan (2017, Section 2.2) give
when defining the effects for arbitrary regimes. Element \\k\\ of a
regime is the exposure at the \\k\\-th distinct value of `time_var`: for
times {1, 3, 4, 6} a regime has four elements. A regime in which the
exposure follows its fitted model outside a window is a dynamic regime,
a different estimand, and is not available. Which pair answers a given
question is the analyst's judgement; only \\a = a^\*\\ is rejected.

**How the mediator is set under `"I"`.** The model is first simulated
once under each regime with the mediator following its fitted model, and
every simulated subject's mediator trajectory \\M(1{:}T)\\ is stored in
a pool. Each decomposition intervention (`Phi00`, `Phi10`, `Phi11`, and
`Phi1_k` with several mediators) then gives each subject the whole
trajectory of a randomly permuted pool member, each mediator permuted
independently, and averages `n_vw` such permutations. The permutation
runs over the whole simulated cohort, so the draw is marginal over the
baseline covariates as well as the confounder history: the
whole-population draw \\G_a\\ of VanderWeele & Tchetgen Tchetgen (2017),
as in the algorithm of Lin et al. (2017, *Epidemiology*). Yamamuro et
al. (2021, Eq. 2) and Lin et al. (2017, *Stat Med*, Eq. 4) are written
for the draw conditional on baseline covariates, but the algorithms of
both papers, and the SAS macros distributed with them, permute the whole
cohort. The two forms agree when the counterfactual mediator
distribution does not depend on the baseline covariates, or when the
outcome mean is additively separable in them and the mediator; otherwise
they are different functionals.

**How the mediator is set under `"N"`.** The mediator model is
re-evaluated on each subject's own simulated history with the exposure
history set to the other regime (Zheng & van der Laan 2017, Eq. 5). Only
the exposure and its first-order lags (`in_recode` entries that copy the
exposure, e.g. `recodes(lag1_A = A)`) are switched, so a mediator model
whose formula or `subset` reads any other exposure-derived column (a
longer lag, a cumulative count, an `out_recode` copy, a column created
by a model's own `recode`) is rejected with an error.

**Total effect and residual.** The total effect is the natural plug-in
contrast `nat1 - nat0`, from runs in which the exposure follows each
regime and the mediator follows its fitted model. Under `"I"` the direct
and indirect effects sum to the interventional overall effect
`Phi11 - Phi00` rather than to the total effect, and the difference is
reported as the decomposition residual. Under `"N"` they sum to the
total effect and no residual is reported.

**Survival outcomes and censoring.** The algorithm is that of Lin et al.
(2017, *Stat Med*, Section 4): models are fitted among those at risk, no
simulated subject is removed, and the cumulative risk \\1-\prod_t(1-\hat
p_t)\\ is accumulated from the predicted hazards. A `"censor"` model
sets the censoring indicator to zero under every intervention, so every
risk is the risk under eliminated loss to follow-up; with
`estimator = "gcomp"` it does not change any reported risk. It marks the
analysis as a survival setting and is used by `estimator = "tmle"`.

**Warnings.** Model-fitting warnings (e.g. non-convergence) are
collected and printed once, with a repeat count, when the function
returns.

## References

Avin, C., Shpitser, I., & Pearl, J. (2005). Identifiability of
path-specific effects. *Proceedings of the 19th International Joint
Conference on Artificial Intelligence*, 357-363.

Lin, S. H., Young, J. G., Logan, R., & VanderWeele, T. J. (2017).
Mediation analysis for a survival outcome with time-varying exposures,
mediators, and confounders. *Statistics in Medicine*, 36(26), 4153–4166.
[doi:10.1002/sim.7426](https://doi.org/10.1002/sim.7426)

Lin, S. H., Young, J., Logan, R., Tchetgen Tchetgen, E. J., &
VanderWeele, T. J. (2017). Parametric mediational g-formula approach to
mediation analysis with time-varying exposures, mediators, and
confounders. *Epidemiology*, 28(2), 266–274.
[doi:10.1097/EDE.0000000000000609](https://doi.org/10.1097/EDE.0000000000000609)

Miles, C. H. (2023). On the causal interpretation of randomised
interventional indirect effects. *Journal of the Royal Statistical
Society: Series B*, 85(4), 1154–1172.
[doi:10.1093/jrsssb/qkad066](https://doi.org/10.1093/jrsssb/qkad066)

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

## See also

[`spec_model`](https://adayim.github.io/causalMed/reference/spec_model.md),
[`recodes`](https://adayim.github.io/causalMed/reference/recodes.md),
[`gformula`](https://adayim.github.io/causalMed/reference/gformula.md)
for total effects, and
[`vignette("causalMed-02-mediation")`](https://adayim.github.io/causalMed/articles/causalMed-02-mediation.md).

## Examples

``` r
data(nonsurvivaldata)

# Models in the order the variables are generated: A -> L1 -> L2 -> M -> Y
models <- list(
  spec_model(A ~ V + lag1_A + lag1_L1 + lag1_L2 + time,
             var_type = "binary", mod_type = "exposure"),
  spec_model(L1 ~ V + A + lag1_L1 + time,
             var_type = "normal", mod_type = "covariate"),
  spec_model(L2 ~ V + A + lag1_L2 + time,
             var_type = "binary", mod_type = "covariate"),
  spec_model(M ~ V + A + L1 + L2 + lag1_M + time,
             var_type = "normal", mod_type = "mediator"),
  spec_model(Y_bin ~ V + A + M + L1 + L2,
             var_type = "binary", mod_type = "outcome")
)

fit <- mediation(
  data = nonsurvivaldata, id_var = "id", time_var = "time",
  base_vars = "V", exposure = "A", outcome = "Y_bin", models = models,
  init_recode = recodes(lag1_A = 0, lag1_L1 = 0, lag1_L2 = 0, lag1_M = 0),
  in_recode   = recodes(lag1_A = A, lag1_L1 = L1, lag1_L2 = L2, lag1_M = M),
  mc_sample = 2000,
  R = 1,          # use R > 1 (e.g. 500) for bootstrap confidence intervals
  quiet = TRUE
)
fit
#> Call:
#> mediation(data = nonsurvivaldata, id_var = "id", base_vars = "V", 
#>     exposure = "A", outcome = "Y_bin", time_var = "time", models = models, 
#>     init_recode = recodes(lag1_A = 0, lag1_L1 = 0, lag1_L2 = 0, 
#>         lag1_M = 0), in_recode = recodes(lag1_A = A, lag1_L1 = L1, 
#>         lag1_L2 = L2, lag1_M = M), mc_sample = 2000, R = 1, quiet = TRUE)
#> 
#> --- Analysis setup ---
#>   Exposure     : A
#>   Mediator(s)  : M
#>   Outcome      : Y_bin  [mean outcome at t = 4, end of follow-up]
#>   Time variable: time  (5 time points: 0 ... 4)
#>   ID variable  : id
#>   Baseline vars: V
#>   Data         : 3,000 individuals, 15,000 observations
#>   Observed subjects following a  (1 1 1 1 1): 1,177 of 3,000 (39.2%)
#>   Observed subjects following a* (0 0 0 0 0): 5 of 3,000 (0.2%)
#>   MC sample    : 2000
#>   Bootstrap R  : none
#>   n_vw         : 2  (permutation draws averaged per pool-drawing intervention)
#>   Seed         : 12345
#>   Mediation    : Interventional effects (IDE/IIE) -- Lin et al. (2017)
#> 
#> --- Marginal mean outcome per intervention --- 
#>   Under interventional effects, each intervention draws its mediators from independently-permuted pools (G):
#>   Phi11 = E[Y(a=1, G1)]:  exposure=1, mediators ~ a=1 pool  [reference]
#>   Phi10 = E[Y(a=1, G0)]:  exposure=1, mediators ~ a=0 pool  [cross-regime]
#>   Phi00 = E[Y(a=0, G0)]:  exposure=0, mediators ~ a=0 pool  [reference]
#>   nat1/nat0 = E[Y(a=1)]/E[Y(a=0)]:  exposure fixed, mediators natural (used for the total effect)
#>    Intervention    Est
#>          <char>  <num>
#> 1:         nat0 0.0934
#> 2:         nat1 0.2433
#> 3:        Phi00 0.0734
#> 4:        Phi10 0.1585
#> 5:        Phi11 0.2207
#>   Observed (nonparametric) mean of Y_bin at t = 4 (end of follow-up): 0.2333
#>   (informal benchmark; interventions fix the exposure, so exact agreement is not expected)
#> 
#> --- Effect decomposition --- 
#>   Direct effect (IDE)   = Phi10 - Phi00
#>   Indirect effect (IIE) = Phi11 - Phi10   (sequential per mediator when N>=2)
#>   IDE + IIE             = Phi11 - Phi00    (interventional overall effect)
#>   Total effect (TE)     = nat1 - nat0      (natural plug-in g-formula)
#>   TE - (Direct+Indirect)= natural TE minus interventional overall effect
#>   Mediation Prop.       = Indirect / (Direct + Indirect)   (percentage)
#>     i.e. a share of the interventional overall effect, NOT of the total
#>     effect; Direct and Mediation Prop. sum to 100%
#>   RD = risk difference;  RR = risk ratio
#>                      Effect      RD     RR
#>                      <char>   <num>  <num>
#> 1:          Indirect effect  0.0622 1.3923
#> 2:            Direct effect  0.0851 2.1599
#> 3:             Total effect  0.1499 2.6047
#> 4: TE - (Direct + Indirect)  0.0026     NA
#> 5:     Mediation Proportion 42.2164     NA
```
