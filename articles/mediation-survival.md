# Mediation Analysis for Survival Outcomes

## Introduction

The mediational g-formula of Lin et al. (2017) was developed for
**survival outcomes** with time-varying exposures, mediators, and
confounders. This vignette walks through that setting:

1.  a single-mediator interventional analysis on a survival outcome,
2.  how to read every row of the output,
3.  the `n_vw` permutation-averaging parameter,
4.  multiple mediators — replicating the Yamamuro et al. (2021)
    simulation setting, with estimates checked against documented true
    values,
5.  censoring, and
6.  restricting models with `subset` (absorbing states).

For an introduction to the package (model specification, recode hooks,
[`gformula()`](https://adayim.github.io/causalMed/reference/gformula.md)
for total effects) see
[`vignette("causalMed-overview")`](https://adayim.github.io/causalMed/articles/causalMed-overview.md).

``` r

library(causalMed)
library(data.table)
#> 
#> Attaching package: 'data.table'
#> The following object is masked from 'package:base':
#> 
#>     %notin%
```

------------------------------------------------------------------------

## Survival Data Requirements

With a survival outcome,
[`mediation()`](https://adayim.github.io/causalMed/reference/mediation.md)
(and
[`gformula()`](https://adayim.github.io/causalMed/reference/gformula.md))
model the **discrete-time hazard**: at each time point, a binary model
with `mod_type = "survival"` estimates the probability of the event
among those still at risk. The simulation accumulates these hazards into
a cumulative incidence, $`1 - \prod_t (1 - h_t)`$, so all reported
quantities — the per-intervention estimates and the effect decomposition
— are **risks of the event by the end of follow-up**.

The data must be in long format with one row per subject per period **at
risk**: once the event occurs, no later rows for that subject may be
present.

``` r

data("survivaldata")
dat <- as.data.table(survivaldata)
head(dat, 8)
#>       id  time         V     A         L     M     Y lag1_A lag1_M    lag1_L
#>    <int> <int>     <num> <int>     <num> <num> <int>  <int>  <num>     <num>
#> 1:     1     0 0.4528886     1 0.6372004     1     1      0      0 0.0000000
#> 2:     2     0 0.6211870     1 0.7814770     1     1      0      0 0.0000000
#> 3:     3     0 0.5418330     1 1.0665595     1     1      0      0 0.0000000
#> 4:     4     0 0.5628616     1 0.8975717     0     0      0      0 0.0000000
#> 5:     4     1 0.5628616     1 1.1214325     1     1      1      0 0.8975717
#> 6:     5     0 0.5246275     0 0.4891121     1     0      0      0 0.0000000
#> 7:     5     1 0.5246275     0 0.6388336     1     0      0      1 0.4891121
#> 8:     5     2 0.5246275     1 1.0670677     0     0      0      1 0.6388336
```

`survivaldata` contains 3 000 subjects followed over five periods
(`time` = 0–4):

| Variable                     | Role                                       |
|------------------------------|--------------------------------------------|
| `id`                         | Subject identifier                         |
| `time`                       | Time index                                 |
| `V`                          | Time-fixed baseline covariate              |
| `A`                          | Time-varying binary exposure               |
| `L`                          | Time-varying continuous confounder         |
| `M`                          | Time-varying binary mediator               |
| `Y`                          | Event indicator (1 = event in this period) |
| `lag1_A`, `lag1_L`, `lag1_M` | Previous-period values                     |

The data-generating ordering within each period is **A → L → M → Y**
(see
[`?survivaldata`](https://adayim.github.io/causalMed/reference/survivaldata.md)):
exposure first, then the confounder (affected by current exposure), then
the mediator (affected by current exposure and confounder), then the
hazard.

------------------------------------------------------------------------

## Single Mediator: Interventional Effects

Models are listed in the temporal order above. Note
`mod_type = "survival"` for the outcome; the formula may include
exposure–mediator interaction terms.

``` r

init_s <- recodes(lag1_A = 0, lag1_L = 0, lag1_M = 0)
in_s   <- recodes(lag1_A = A, lag1_L = L, lag1_M = M)

models_surv <- list(
  spec_model(A ~ V + lag1_A + lag1_L + time,
             var_type = "binary", mod_type = "exposure"),
  spec_model(L ~ V + A + lag1_L + time,
             var_type = "normal", mod_type = "covariate"),
  spec_model(M ~ V + A + L + lag1_M + time,
             var_type = "binary", mod_type = "mediator"),
  spec_model(Y ~ V + A + M + L + A:M + time,
             var_type = "binary", mod_type = "survival")   # discrete-time hazard
)
```

``` r

fit_surv <- mediation(
  data           = dat,
  id_var         = "id",
  time_var       = "time",
  base_vars      = "V",
  exposure       = "A",
  outcome        = "Y",
  models         = models_surv,
  init_recode    = init_s,
  in_recode      = in_s,
  mediation_type = "I",
  mc_sample      = 20000,
  R              = 1,      # set R > 1 for bootstrap CIs (see below)
  quiet          = TRUE,
  seed           = 2026
)
```

### Reading the per-intervention table

``` r

fit_surv$effect_size
#>    Intervention       Est
#>          <char>     <num>
#> 1:         nat0 0.8369068
#> 2:         nat1 0.9753288
#> 3:        Phi00 0.8369122
#> 4:        Phi10 0.9623516
#> 5:        Phi11 0.9753319
```

Each row is the simulated cumulative incidence under one intervention:

- **`nat0` / `nat1`** — exposure fixed to 0 / 1, mediator following its
  fitted model given each individual’s own history. Their contrast is
  the natural plug-in **total effect**.
- **`Phi00` / `Phi11`** — exposure fixed to 0 / 1, mediator drawn from
  the *permuted marginal pool* collected under that same exposure level
  ($`E[Y_{0,G_0}]`$, $`E[Y_{1,G_1}]`$). These are the interventional
  references.
- **`Phi10`** — the cross-world intervention: exposure fixed to 1,
  mediator drawn from the a = 0 pool ($`E[Y_{1,G_0}]`$).

### Reading the decomposition

``` r

fit_surv$estimate
#>                                   Effect           RD       RR
#>                                   <char>        <num>    <num>
#> 1:                       Indirect effect 1.298035e-02 1.013488
#> 2:                         Direct effect 1.254393e-01 1.149884
#> 3:                          Total effect 1.384220e-01 1.165397
#> 4:              TE - (Direct + Indirect) 2.329146e-06       NA
#> 5:                  Mediation Proportion 9.379053e+00       NA
#> 6: Mediation Proportion (multiplicative) 9.377529e+00       NA
```

Row by row:

- **Indirect effect (IIE)** = `Phi11 − Phi10`: the change in risk from
  shifting the *population distribution* of the mediator from its
  never-treated to its always-treated form, while exposure is held at 1.
- **Direct effect (IDE)** = `Phi10 − Phi00`: the effect of exposure with
  the mediator distribution held at its never-treated form.
- **Total effect (TE)** = `nat1 − nat0`: the ordinary g-formula total
  effect.
- **TE − (Direct + Indirect)**: IDE + IIE sum to the *interventional
  overall effect* `Phi11 − Phi00`, not to TE; this row is the difference
  (a mediated-interaction residual). It is usually small; a large value
  signals strong exposure–mediator interaction in the outcome process.
  (For `mediation_type = "N"` the decomposition is exact and this row is
  absent.)
- **Mediation Proportion** = (TE − IDE) / TE × 100: the share of the
  total effect *not* acting directly, i.e. the indirect effects plus the
  residual (Yamamuro et al. 2021).
- **Mediation Proportion (multiplicative)** =
  $`RR_{IDE}(RR_{IIE}-1)/(RR_{OE}-1) \times 100`$ on the risk-ratio
  scale (Lin et al. 2017, Table 2).

The `RD` column is the risk-difference scale, `RR` the risk-ratio scale
(`RR` is not applicable to the residual and proportion rows). With
`R > 1` both scales get bootstrap standard errors and percentile/normal
confidence intervals.

`print(fit_surv)` displays both tables together with the intervention
definitions as a legend, an analysis-setup summary (mediator, outcome,
data dimensions, `n_vw`), and an observed (nonparametric)
cumulative-incidence benchmark computed directly from the data.

### The `n_vw` argument

Every cross-world intervention draws the mediator by randomly permuting
the pool of simulated mediator trajectories. `n_vw` controls how many
independent permutations are averaged per intervention (default `2`,
matching the SAS `mGFORMULA` macro). Averaging reduces Monte Carlo noise
from the permutation step at the cost of one extra simulation pass per
intervention; `n_vw = 1` is faster and slightly noisier. It has no
effect on natural effects (`mediation_type = "N"`), which do not use
permutation.

``` r

mediation(..., n_vw = 1)   # single permutation per intervention: faster, noisier
```

------------------------------------------------------------------------

## Multiple Mediators: the Yamamuro et al. (2021) Simulation

With `mediation_type = "I"`, any number of mediators can be analysed by
supplying several `mod_type = "mediator"` models **in temporal order**.
We demonstrate this on `yamamurodata`, a dataset simulated from the
data-generating process of the Yamamuro et al. (2021) simulation study —
a time-varying treatment `A`, a confounder `L`, two sequential mediators
`M1` and `M2`, and a survival outcome over three visits, ordered **A → L
→ M1 → M2 → Y** within each visit. Because the process is known, the
**true interventional effects are documented in
[`?yamamurodata`](https://adayim.github.io/causalMed/reference/yamamurodata.md)**
(computed at $`n = 10^7`$), so the estimates below can be checked
against the truth.

``` r

data("yamamurodata")
yam <- as.data.table(yamamurodata)
head(yam, 6)
#>       id  time     V     A        L       M1       M2     Y lag1_A   lag1_L
#>    <int> <int> <int> <int>    <num>    <num>    <num> <int>  <int>    <num>
#> 1:     1     0     0     0 22.99932 154.1412 6560.201     0      0  0.00000
#> 2:     1     1     0     0 22.72770 150.8244 6674.132     0      0 22.99932
#> 3:     1     2     0     1 23.66563 142.4219 6724.582     0      0 22.72770
#> 4:     2     0     1     0 24.95263 188.5737 6593.398     0      0  0.00000
#> 5:     2     1     1     1 24.49130 141.9657 6196.278     0      0 24.95263
#> 6:     2     2     1     1 24.17914 132.8109 5922.934     0      1 24.49130
#>     lag1_M1  lag1_M2   L0base  M10base  M20base
#>       <num>    <num>    <num>    <num>    <num>
#> 1:   0.0000    0.000 22.99932 154.1412 6560.201
#> 2: 154.1412 6560.201 22.99932 154.1412 6560.201
#> 3: 150.8244 6674.132 22.99932 154.1412 6560.201
#> 4:   0.0000    0.000 24.95263 188.5737 6593.398
#> 5: 188.5737 6593.398 24.95263 188.5737 6593.398
#> 6: 141.9657 6196.278 24.95263 188.5737 6593.398
```

Three specification details are worth noting:

- **Visit indicators.** The published models give each visit its own
  intercept. Write these as `I(as.integer(time == k))` rather than
  `factor(time)`: during simulation each Monte Carlo time slice carries
  a single `time` value, so a factor would drop the unobserved levels
  while the numeric indicator predicts safely at every step.
- **`subset = time > 0`.** The visit-0 values of `A`, `L`, `M1`, `M2`
  are baseline draws, not model output, so the time-varying models are
  fitted and simulated only for `time > 0`. Visit 0 is seeded in the
  Monte Carlo cohort from the time-fixed baseline columns (`L0base`,
  `M10base`, `M20base`) via `init_recode`.
- The models include the quadratic and lag terms of the published
  “correctly specified” scenario.

``` r

models_yam <- list(
  spec_model(A ~ V + lag1_A + lag1_L + I(lag1_L^2) + lag1_M1 + lag1_M2 +
               I(as.integer(time == 2)),
             var_type = "binary", mod_type = "exposure", subset = time > 0),
  spec_model(L ~ V + A + lag1_L + I(lag1_L^2) + lag1_M1 + lag1_M2 +
               I(as.integer(time == 2)),
             var_type = "normal", mod_type = "covariate", subset = time > 0),
  spec_model(M1 ~ V + A + L + I(L^2) + lag1_L + I(lag1_L^2) + lag1_M1 + lag1_M2 +
               I(as.integer(time == 2)),
             var_type = "normal", mod_type = "mediator", subset = time > 0),
  spec_model(M2 ~ V + A + L + I(L^2) + lag1_L + I(lag1_L^2) + lag1_M1 + lag1_M2 +
               I(as.integer(time == 2)),
             var_type = "normal", mod_type = "mediator", subset = time > 0),
  spec_model(Y ~ V + A + L + I(L^2) + M1 + M2 + M1:M2 +
               I(as.integer(time == 1)) + I(as.integer(time == 2)),
             var_type = "binary", mod_type = "survival")
)

fit_yam <- mediation(
  data           = yam,
  id_var         = "id",
  time_var       = "time",
  base_vars      = c("V", "L0base", "M10base", "M20base"),
  exposure       = "A",
  outcome        = "Y",
  models         = models_yam,
  init_recode    = recodes(L = L0base, M1 = M10base, M2 = M20base,
                           lag1_A = 0, lag1_L = 0, lag1_M1 = 0, lag1_M2 = 0),
  in_recode      = recodes(lag1_A = A, lag1_L = L, lag1_M1 = M1, lag1_M2 = M2),
  mediation_type = "I",
  n_vw           = 2,      # matches the SAS mGFORMULA macro
  mc_sample      = 20000,
  R              = 1,
  quiet          = TRUE,
  seed           = 2026
)

fit_yam$effect_size
#>    Intervention        Est
#>          <char>      <num>
#> 1:         nat0 0.11248807
#> 2:         nat1 0.03509880
#> 3:        Phi00 0.11446276
#> 4:        Phi10 0.06851088
#> 5:       Phi1_1 0.04526965
#> 6:        Phi11 0.03573016
```

For $`N`$ mediators the intervention list grows to $`4 + N`$: the
intermediate interventions `Phi1_k` switch the first $`k`$ mediators to
the a = 1 pool while the rest stay at the a = 0 pool. Each mediator’s
indirect effect is the **sequential contrast**
$`IIE(M_k) = \Phi_{1,k} - \Phi_{1,k-1}`$, labelled
`Indirect effect (<name>)` in the output, and $`IDE + \sum_k IIE(M_k)`$
equals the interventional overall effect `Phi11 − Phi00`. Each
mediator’s pool is permuted independently.

### Comparing against the true values

``` r

truth <- data.table(
  Effect = c("Total effect", "Direct effect", "Indirect effect (M1)",
             "Indirect effect (M2)", "TE - (Direct + Indirect)"),
  True   = c(-6.36, -3.20, -2.29, -0.97, 0.10)   # from ?yamamurodata
)
est <- as.data.table(fit_yam$estimate)[, .(Effect, Estimate = RD * 100)]
merge(truth, est, by = "Effect", sort = FALSE)
#>                      Effect  True   Estimate
#>                      <char> <num>      <num>
#> 1:             Total effect -6.36 -7.7389265
#> 2:            Direct effect -3.20 -4.5951879
#> 3:     Indirect effect (M1) -2.29 -2.3241236
#> 4:     Indirect effect (M2) -0.97 -0.9539488
#> 5: TE - (Direct + Indirect)  0.10  0.1343338
```

The estimates reproduce the structure of the truth: all signs and the
relative magnitudes (IDE \> IIE via M1 \> IIE via M2, small positive
residual) are recovered, and the two indirect effects are estimated
closely. The total and direct effects deviate by roughly 1–1.5
percentage points — sampling error from fitting on a single dataset. The
published study averages 1000 replicate datasets, under which the mean
estimates match the true values (Yamamuro et al. 2021, and reproduced
with this package during development); a per-dataset bootstrap (`R > 1`)
quantifies this uncertainty in applied use.

Because the decomposition is sequential, the *order* of the mediator
models matters: it must reflect the assumed temporal/causal ordering
among the mediators. Multiple mediators are not available for
`mediation_type = "N"` (the natural-effects references define a single
mediator only);
[`mediation()`](https://adayim.github.io/causalMed/reference/mediation.md)
stops with an error in that case.

------------------------------------------------------------------------

## Censoring

If follow-up can end for reasons other than the event, include a
censoring indicator with `mod_type = "censor"`. To illustrate, we censor
some follow-up in `survivaldata` with probability depending on the
confounder (informative censoring), ending each subject’s follow-up at
their first censoring:

``` r

set.seed(11)
dat_c <- copy(dat)
dat_c[, C := rbinom(.N, 1, plogis(-3.5 + 0.8 * L))]
dat_c[, after_cens := cumsum(shift(C, fill = 0)), by = id]
dat_c <- dat_c[after_cens == 0][, after_cens := NULL]
dat_c[C == 1, Y := 0L]   # censored before the event in that period
```

Add a censoring model to the list (after the covariates/mediator it
depends on, before the survival model):

``` r

models_cens <- list(
  spec_model(A ~ V + lag1_A + lag1_L + time,
             var_type = "binary", mod_type = "exposure"),
  spec_model(L ~ V + A + lag1_L + time,
             var_type = "normal", mod_type = "covariate"),
  spec_model(M ~ V + A + L + lag1_M + time,
             var_type = "binary", mod_type = "mediator"),
  spec_model(C ~ V + A + L + time,
             var_type = "binary", mod_type = "censor"),     # censoring model
  spec_model(Y ~ V + A + M + L + A:M + time,
             var_type = "binary", mod_type = "survival")
)

fit_cens <- mediation(
  data           = dat_c,
  id_var         = "id",
  time_var       = "time",
  base_vars      = "V",
  exposure       = "A",
  outcome        = "Y",
  models         = models_cens,
  init_recode    = init_s,
  in_recode      = in_s,
  mediation_type = "I",
  mc_sample      = 20000,
  R              = 1,
  quiet          = TRUE,
  seed           = 2026
)

fit_cens$estimate
#>                                   Effect            RD       RR
#>                                   <char>         <num>    <num>
#> 1:                       Indirect effect  1.496432e-02 1.015768
#> 2:                         Direct effect  1.298705e-01 1.158545
#> 3:                          Total effect  1.448257e-01 1.176800
#> 4:              TE - (Direct + Indirect) -9.101517e-06       NA
#> 5:                  Mediation Proportion  1.032636e+01       NA
#> 6: Mediation Proportion (multiplicative)  1.033199e+01       NA
```

In the simulation, censoring is **abolished** in every intervention (all
interventions fix the exposure, and the censoring indicator is set to
0), so the reported risks are counterfactual risks *in the absence of
censoring* — the standard g-formula treatment of right-censoring. The
censoring model’s role is to let the hazard model be fitted on data
where censoring depends on measured covariates. Here the estimates from
the censored data are close to the full-data estimates above, as
expected when censoring depends only on modelled covariates.

------------------------------------------------------------------------

## Restricting Models with `subset` (Absorbing States)

`spec_model(subset = ...)` fits and simulates a model only on rows
meeting a condition. The classic use is an **absorbing state** — for
example, in the GvHD analysis of Keil et al. (2014) the exposure can
only switch on once, so its model is estimated among the not-yet-exposed
and the value is carried forward deterministically afterwards:

``` r

# exposure can occur only while gvhdm1 == 0 …
spec_model(gvhd ~ all + cmv + male + age + ...,
           var_type = "binary", mod_type = "exposure",
           subset = gvhdm1 == 0)

# … and is locked at 1 afterwards via the end-of-step hook
out_recode = recodes(gvhd = ifelse(gvhdm1 == 1, 1, gvhd))
```

Rows excluded by `subset` keep their current value at that step, so the
`out_recode` carry-forward completes the absorbing behaviour. The
package ships the `gvhd` dataset (see
[`?gvhd`](https://adayim.github.io/causalMed/reference/gvhd.md)) used in
that paper, and a complete total-effect g-formula analysis on it — three
absorbing states, a censoring model, and restricted cubic splines — is
given in `references/paper_example.R` and walked through in
[`vignette("causalMed-overview")`](https://adayim.github.io/causalMed/articles/causalMed-overview.md).

------------------------------------------------------------------------

## Bootstrap Confidence Intervals

As elsewhere in the package, set `R > 1` for subject-level bootstrap CIs
and optionally register a parallel plan first:

``` r

library(future)
plan(multisession)

fit_ci <- mediation(
  data           = dat,
  id_var         = "id",
  time_var       = "time",
  base_vars      = "V",
  exposure       = "A",
  outcome        = "Y",
  models         = models_surv,
  init_recode    = init_s,
  in_recode      = in_s,
  mediation_type = "I",
  mc_sample      = 20000,
  R              = 500,
  seed           = 2026
)

plan(sequential)

# estimate now carries Sd, percentile and normal CIs on the RD and RR scales
fit_ci$estimate

# the per-replicate bootstrap draws are also retained, for custom diagnostics
# (e.g. counting non-finite Mediation Proportion replicates):
fit_ci$boot_estimates$effects
```

------------------------------------------------------------------------

## Natural Effects and Intermediate Confounding

`mediation_type = "N"` (Zheng & van der Laan 2017) is also defined for
survival outcomes. However, in this dataset the confounder `L` is
affected by the current exposure — an *intermediate confounder* — and
natural direct and indirect effects are **not identifiable** in that
setting.
[`mediation()`](https://adayim.github.io/causalMed/reference/mediation.md)
detects the exposure on the right-hand side of a covariate model and
warns:

``` r

fit_nat <- mediation(
  data           = dat,
  id_var         = "id",
  time_var       = "time",
  base_vars      = "V",
  exposure       = "A",
  outcome        = "Y",
  models         = models_surv,
  init_recode    = init_s,
  in_recode      = in_s,
  mediation_type = "N",
  mc_sample      = 20000,
  R              = 1,
  quiet          = TRUE,
  seed           = 2026
)
#> Warning: mediation_type = "N" requested, but covariate model(s) for {L} include
#> the exposure 'A' on the right-hand side, indicating an intermediate
#> (exposure-affected) confounder. Natural direct and indirect effects are NOT
#> identifiable from observational data in this setting (Avin, Shpitser & Pearl
#> 2005; VanderWeele 2014; VanderWeele & Tchetgen Tchetgen 2017). Consider
#> mediation_type = "I" for identifiable interventional (randomized-analogue)
#> effects.
```

The interventional effects (`"I"`) remain identifiable under
intermediate confounding and are the recommended estimand here — which
is why they are the package default.

------------------------------------------------------------------------

## References

- Lin, S.-H., Young, J. G., Logan, R., & VanderWeele, T. J. (2017).
  Mediation analysis for a survival outcome with time-varying exposures,
  mediators, and confounders. *Statistics in Medicine*, 36, 4153–4166.
- VanderWeele, T. J., & Tchetgen Tchetgen, E. J. (2017). Mediation
  analysis with time varying exposures and mediators. *Journal of the
  Royal Statistical Society: Series B*, 79(3), 917–938.
- Yamamuro, S., Shinozaki, T., Iimuro, S., & Matsuyama, Y. (2021).
  Mediational g-formula for time-varying treatment and repeated-measured
  multiple mediators. *Statistical Methods in Medical Research*, 30(8),
  1782–1799.
- Zheng, W., & van der Laan, M. (2017). Longitudinal mediation analysis
  with time-varying mediators and exposures, with application to
  survival outcomes. *Journal of Causal Inference*, 5(2).
- Keil, A. P., Edwards, J. K., Richardson, D. B., Naimi, A. I., &
  Cole, S. R. (2014). The parametric g-formula for time-to-event data:
  intuition and a worked example. *Epidemiology*, 25(6), 889–897.
- Westreich, D., Cole, S. R., Young, J. G., et al. (2012). The
  parametric g-formula to estimate the effect of highly active
  antiretroviral therapy on incident AIDS or death. *Statistics in
  Medicine*, 31, 2000–2009.
