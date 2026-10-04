# Causal Mediation Analysis with causalMed

## Introduction

[`mediation()`](https://adayim.github.io/causalMed/reference/mediation.md)
splits the effect of a time-varying exposure into a **direct** effect
and an **indirect** effect through one or more time-varying mediators,
including when confounders are affected by earlier exposure, the setting
the mediational g-formula was developed for (Lin et al. 2017;
VanderWeele & Tchetgen Tchetgen 2017). This vignette works through a
survival example and reads every row of the output, then covers other
exposure regimes, multiple mediators, censoring, bootstrap intervals and
natural effects.

Data format,
[`spec_model()`](https://adayim.github.io/causalMed/reference/spec_model.md)
and the
[`recodes()`](https://adayim.github.io/causalMed/reference/recodes.md)
hooks are introduced in
[`vignette("causalMed-01-overview")`](https://adayim.github.io/causalMed/articles/causalMed-01-overview.md).

``` r

library(causalMed)
library(data.table)
```

## The two estimands

`mediation_type` chooses between two definitions of the direct and
indirect effects.

**Natural effects** (`"N"`; Zheng & van der Laan 2017) decompose the
total effect exactly:

``` math
\underbrace{E[Y_{1,M(1)}] - E[Y_{0,M(0)}]}_{\text{Total effect}} =
  \underbrace{E[Y_{1,M(0)}] - E[Y_{0,M(0)}]}_{\text{Direct effect}} +
  \underbrace{E[Y_{1,M(1)}] - E[Y_{1,M(0)}]}_{\text{Indirect effect}}
```

Here the mediator is drawn from its **conditional** distribution: at
each time point it is predicted from the subject’s own history, with the
exposure and its first-order lags set to the other level.

**Interventional effects** (`"I"`, the default; Lin et al. 2017) set the
mediator to a random draw $`G_{a^*}`$ from its **marginal** distribution
under $`a^*`$, drawn independently for each mediator. The direct and
indirect effects then sum to the **interventional overall effect**
$`E[Y_{1,G_1}] - E[Y_{0,G_0}]`$, which generally differs from the total
effect $`E[Y_1] - E[Y_0]`$.
[`mediation()`](https://adayim.github.io/causalMed/reference/mediation.md)
reports the total effect separately and adds the difference as a row,
`TE - (Direct + Indirect)`, the decomposition residual.

|  | Interventional (`"I"`) | Natural (`"N"`) |
|----|----|----|
| Mediator drawn from | Marginal distribution (permutation) | Conditional distribution (own history) |
| Cross-world assumption needed | No | Only to read the effects as individual-level |
| Identified with exposure-affected confounders | Yes | Yes, as the conditional-draw estimand |
| Several mediators | Yes (Yamamuro et al. 2021) | No |
| Reference | Lin et al. (2017) | Zheng & van der Laan (2017) |

Zheng & van der Laan (2017, Lemma 1) identify the `"N"` effects under
sequential randomization and positivity. Reading them as
**individual-level** natural effects, contrasts of each subject’s own
counterfactual mediator, additionally requires a cross-world assumption
that is not expected to hold when a mediator-outcome confounder is
itself affected by the exposure (Avin, Shpitser & Pearl 2005;
VanderWeele & Tchetgen Tchetgen 2017). For that setting VanderWeele &
Tchetgen Tchetgen (2017) propose the interventional effects, which is
why `"I"` is the default. With `"N"`,
[`mediation()`](https://adayim.github.io/causalMed/reference/mediation.md)
warns when a covariate model includes the exposure; the check reads
formulas only and cannot confirm or rule out such a confounder.

## Requirements

- **Binary exposure**, coded 0/1.
- **At least one mediator model** (`mod_type = "mediator"`). Several
  mediators, in temporal order, are supported under `"I"` only.
- **Model order** follows the order in which variables are generated
  within a time point, typically `A → L → M → Y`, or `A → M → L → Y`
  when the confounders respond to the mediator (the setting of Lin et
  al. 2017). A warning is given if the exposure follows the mediator or
  the mediator follows the outcome.
- **Outcome**: binary or continuous at the end of follow-up
  (`mod_type = "outcome"`) or a discrete-time event
  (`mod_type = "survival"`).

## Survival data

With a survival outcome, a `mod_type = "survival"` model estimates the
discrete-time hazard among subjects still at risk. The simulation
accumulates it into the cumulative incidence $`1 - \prod_t (1 - h_t)`$,
so every reported quantity is a **risk of the event by the end of
follow-up**. The data must have one row per subject per period at risk,
with no rows after the event.

``` r

data("survivaldata")
dat <- as.data.table(survivaldata)
head(dat, 8)
#>       id  time          V     A           L     M     Y     C lag1_A lag1_M
#>    <int> <int>      <num> <int>       <num> <num> <int> <int>  <int>  <num>
#> 1:     1     0 -0.9422279     1  0.10013047     0     0     0      0      0
#> 2:     1     1 -0.9422279     0 -0.05257470     0     0     0      1      0
#> 3:     1     2 -0.9422279     0  0.02384068     1     0     0      0      0
#> 4:     1     3 -0.9422279     0  0.14223027     0     0     0      0      1
#> 5:     1     4 -0.9422279     0  0.30567376     0     0     0      0      0
#> 6:     2     0  2.4237402     1  2.16300892     1     0     0      0      0
#> 7:     2     1  2.4237402     1  2.63056349     0     0     0      1      1
#> 8:     2     2  2.4237402     1  2.95539778     1     0     0      1      0
#>         lag1_L
#>          <num>
#> 1:  0.00000000
#> 2:  0.10013047
#> 3: -0.05257470
#> 4:  0.02384068
#> 5:  0.14223027
#> 6:  0.00000000
#> 7:  2.16300892
#> 8:  2.63056349
```

`survivaldata` follows 3,000 subjects over five periods (`time` 0 to 4):

| Variable                     | Role                   |
|------------------------------|------------------------|
| `id`, `time`                 | Subject and time index |
| `V`                          | Baseline covariate     |
| `A`                          | Binary exposure        |
| `L`                          | Continuous confounder  |
| `M`                          | Binary mediator        |
| `Y`                          | Event indicator        |
| `C`                          | Loss to follow-up      |
| `lag1_A`, `lag1_L`, `lag1_M` | Previous-period values |

Within each period the order is **A → L → M → Y** (see
[`?survivaldata`](https://adayim.github.io/causalMed/reference/survivaldata.md)).

## Single mediator: interventional effects

The models follow that order, with a censoring model for loss to
follow-up (see *Censoring* below) and a hazard model that may include
exposure-mediator interactions:

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
  spec_model(C ~ V + L + time,
             var_type = "binary", mod_type = "censor"),    # loss to follow-up
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

`print(fit_surv)` shows everything below together, with a legend of the
interventions, the analysis setup, and an observed (nonparametric)
cumulative-incidence benchmark.

### The interventions

``` r

fit_surv$effect_size
#>    Intervention       Est
#>          <char>     <num>
#> 1:         nat0 0.4106340
#> 2:         nat1 0.7868999
#> 3:        Phi00 0.4177377
#> 4:        Phi10 0.5903442
#> 5:        Phi11 0.7823181
```

Each row is the simulated cumulative incidence under one intervention:

- **`nat0` / `nat1`**: exposure fixed to the reference regime $`a^*`$ /
  the exposure regime $`a`$ (by default never / always exposed),
  mediator following its fitted model. Their contrast is the **total
  effect**.
- **`Phi00` / `Phi11`**: exposure fixed to $`a^*`$ / $`a`$, mediator
  drawn from the permuted pool of mediator trajectories simulated under
  that same regime.
- **`Phi10`**: exposure fixed to $`a`$, mediator drawn from the $`a^*`$
  pool.

The pools hold each simulated subject’s whole mediator trajectory; each
subject is given the trajectory of a randomly permuted pool member (see
[`?mediation`](https://adayim.github.io/causalMed/reference/mediation.md)).

### The decomposition

``` r

fit_surv$estimate
#>                      Effect          RD       RR
#>                      <char>       <num>    <num>
#> 1:          Indirect effect  0.19197390 1.325190
#> 2:            Direct effect  0.17260655 1.413194
#> 3:             Total effect  0.37626584 1.916305
#> 4: TE - (Direct + Indirect)  0.01168539       NA
#> 5:     Mediation Proportion 52.65611453       NA
```

- **Indirect effect** = `Phi11 − Phi10`: with exposure, the change in
  risk from shifting the mediator’s distribution from its unexposed to
  its exposed form.
- **Direct effect** = `Phi10 − Phi00`: the effect of exposure with the
  mediator drawn from its unexposed distribution.
- **Total effect** = `nat1 − nat0`: the ordinary g-formula total effect.
- **TE − (Direct + Indirect)**: the direct and indirect effects sum to
  the interventional overall effect `Phi11 − Phi00`, not to the total
  effect; this row is the difference. It is absent under `"N"`.
- **Mediation Proportion**: the indirect effect as a percentage of the
  interventional overall effect, the quantity Lin et al. (2017, Table 2)
  report. It is not a share of the total effect.

`RD` is the risk difference and `RR` the risk ratio (`NA` for the last
two rows). With `R > 1` both gain bootstrap standard errors and
confidence limits.

### Permutation averaging: `n_vw`

Each pool-drawing intervention (`Phi00`, `Phi10`, `Phi11`) averages
`n_vw` independent permutations of the pool (default 2, as in the SAS
`mGFORMULA` macro). `n_vw = 1` is faster, with more Monte Carlo noise.
It has no effect under `"N"`, which does not permute.

``` r

mediation(..., n_vw = 1)
```

### Other exposure regimes

By default the effects compare always exposed with never exposed, the
example Lin et al. (2017, Section 2.2) and Zheng & van der Laan (2017,
Section 2.2) give when defining the effects for arbitrary regimes.
`exposure_regime` and `reference_regime` take any two static regimes as
0/1 vectors, one element per distinct time point (a scalar is recycled).
Exposure at times 1 and 2 only, against never exposed:

``` r

fit_window <- mediation(
  data             = dat,
  id_var           = "id",
  time_var         = "time",
  base_vars        = "V",
  exposure         = "A",
  outcome          = "Y",
  models           = models_surv,
  init_recode      = init_s,
  in_recode        = in_s,
  mediation_type   = "I",
  exposure_regime  = c(0, 1, 1, 0, 0),   # exposed at t = 1 and 2 only
  reference_regime = 0,                  # never exposed
  mc_sample        = 20000,
  R                = 1,
  quiet            = TRUE,
  seed             = 2026
)
fit_window$estimate
#>                      Effect           RD       RR
#>                      <char>        <num>    <num>
#> 1:          Indirect effect  0.111479530 1.226072
#> 2:            Direct effect  0.075377855 1.180443
#> 3:             Total effect  0.192451782 1.468670
#> 4: TE - (Direct + Indirect)  0.005594397       NA
#> 5:     Mediation Proportion 59.660221446       NA
```

The rows read as before, with $`a`$ now the windowed regime. Two points:

- The pools hold whole mediator trajectories, so the indirect effect
  includes the window’s effect on the mediator at $`t`$ = 3 and 4, after
  the window has closed.
- The zeros outside the window force the exposure to 0 there; they do
  not let it follow its fitted model. That would be a dynamic regime, a
  different estimand, and is not available.

Which pair of regimes answers a given question is the analyst’s
judgement; the package checks only that they are 0/1, of the right
length, and different. [`print()`](https://rdrr.io/r/base/print.html)
also reports how many observed subjects follow each regime:

``` r

fit_window$data_summary$regime_support
#>   regime    values n_following prop_following n_complete
#> 1      a 0 1 1 0 0         489          0.163         23
#> 2     a* 0 0 0 0 0         528          0.176        162
```

`n_following` counts subjects whose exposure matches the regime at every
time they were observed; `n_complete` counts those among subjects
observed at every time point. The estimation does not use these counts:
the g-formula produces an estimate whether or not anyone followed the
regime.

## Multiple mediators: the Yamamuro et al. (2021) simulation

Under `"I"`, several mediators are analysed by listing their models in
temporal order. `yamamurodata` was simulated from the data-generating
process of Yamamuro et al. (2021): a treatment `A`, a confounder `L`,
two sequential mediators `M1` and `M2`, and a survival outcome over
three visits, ordered **A → L → M1 → M2 → Y**. Its true effects are
given in
[`?yamamurodata`](https://adayim.github.io/causalMed/reference/yamamurodata.md).

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

Three details of the specification:

- **Visit indicators** are written `I(as.integer(time == k))` rather
  than `factor(time)`: each simulated time step has a single `time`
  value, so a factor would lose its other levels.
- **`subset = time > 0`**: visit-0 values are baseline draws, so the
  time-varying models apply from visit 1. Visit 0 is set from the
  baseline columns `L0base`, `M10base`, `M20base` in `init_recode`.
- The formulas include the quadratic and lag terms of the published
  correctly specified scenario.

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

With $`N`$ mediators there are $`4 + N`$ interventions: `Phi1_k` draws
the first $`k`$ mediators from the exposed pool and the rest from the
unexposed pool. The indirect effect through mediator $`k`$ is the
**sequential contrast** $`\Phi_{1,k} - \Phi_{1,k-1}`$, labelled
`Indirect effect (<name>)`, and the direct effect plus all indirect
effects equals `Phi11 − Phi00`. Because the decomposition is sequential,
the order of the mediator models matters.

Comparing with the true values (in percentage points):

``` r

truth <- data.table(
  Effect = c("Total effect", "Direct effect", "Indirect effect (M1)",
             "Indirect effect (M2)", "TE - (Direct + Indirect)"),
  True   = c(-6.36, -3.20, -2.29, -0.97, 0.10)   # from ?yamamurodata
)
est <- as.data.table(fit_yam$estimate)[, .(Effect, Estimate = RD * 100)]
cmp <- merge(truth, est, by = "Effect", sort = FALSE)
cmp[, Difference := Estimate - True]
cmp
#>                      Effect  True   Estimate  Difference
#>                      <char> <num>      <num>       <num>
#> 1:             Total effect -6.36 -7.7389265 -1.37892647
#> 2:            Direct effect -3.20 -4.5951879 -1.39518787
#> 3:     Indirect effect (M1) -2.29 -2.3241236 -0.03412365
#> 4:     Indirect effect (M2) -0.97 -0.9539488  0.01605121
#> 5: TE - (Direct + Indirect)  0.10  0.1343338  0.03433385
```

The true values are large-sample values from the data-generating
process. The estimate comes from one dataset of 10,000 subjects and one
Monte Carlo run, so the differences combine sampling and Monte Carlo
error, neither of which this single run quantifies; a bootstrap
(`R > 1`) does.

## Censoring

When follow-up can end before the event, add a censoring indicator with
`mod_type = "censor"`. To illustrate censoring that depends on the
confounder, we add extra loss to follow-up to `survivaldata`, ending
each subject’s follow-up at the first censoring:

``` r

set.seed(11)
dat_c <- copy(dat)
dat_c[, C := rbinom(.N, 1, plogis(-3.5 + 0.8 * L))]
dat_c[, after_cens := cumsum(shift(C, fill = 0)), by = id]
dat_c <- dat_c[after_cens == 0][, after_cens := NULL]
dat_c[C == 1, Y := 0L]   # censored before the event in that period
```

The censoring model goes after the variables it depends on and before
the survival model:

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
#>                      Effect          RD       RR
#>                      <char>       <num>    <num>
#> 1:          Indirect effect  0.19200968 1.344445
#> 2:            Direct effect  0.16502821 1.420542
#> 3:             Total effect  0.37218202 1.974708
#> 4: TE - (Direct + Indirect)  0.01514413       NA
#> 5:     Mediation Proportion 53.77851536       NA
```

Under every intervention the censoring indicator is set to 0, so the
reported risks are risks in the absence of censoring, as in the
g-formula treatment of right-censoring (Robins 1986; Westreich et
al. 2012). With the default `estimator = "gcomp"` the censoring model
does not change these risks: the hazard model is fitted on the rows
still at risk, which identifies them when censoring is independent of
the event given the modelled history, an assumption the package cannot
check. The censoring model is used by the targeted estimator (see
below).

## Absorbing states with `subset`

`spec_model(subset = ...)` fits and simulates a model only on the rows
that meet a condition. The classic use is an **absorbing state**, such
as the GvHD exposure of Keil et al. (2014), which can switch on only
once: its model is fitted among the not-yet-exposed, and `out_recode`
carries the value forward.

``` r

# the exposure can occur only while gvhdm1 == 0 ...
spec_model(gvhd ~ all + cmv + male + age + ...,
           var_type = "binary", mod_type = "exposure",
           subset = gvhdm1 == 0)

# ... and stays at 1 afterwards
out_recode = recodes(gvhd = ifelse(gvhdm1 == 1, 1, gvhd))
```

Rows excluded by `subset` keep their current value. The full GvHD
analysis
([`?gvhd`](https://adayim.github.io/causalMed/reference/gvhd.md)) is in
[`vignette("causalMed-03-gformula")`](https://adayim.github.io/causalMed/articles/causalMed-03-gformula.md).

## Bootstrap confidence intervals

Set `R > 1` for subject-level bootstrap intervals, optionally with a
parallel plan:

``` r

future::plan(future::multisession)

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

future::plan(future::sequential)

fit_ci$estimate                  # with Sd, percentile and normal limits
fit_ci$boot_estimates$effects    # the per-replicate estimates
```

## Natural effects

`mediation_type = "N"` gives the natural effects of Zheng & van der Laan
(2017), for survival outcomes too. Here the confounder `L` responds to
the current exposure, so
[`mediation()`](https://adayim.github.io/causalMed/reference/mediation.md)
warns: the effects are still identified (their Lemma 1), but not as
individual-level natural effects.

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
#> Warning: mediation_type = "N" requested, and covariate model(s) for {L} include
#> the exposure 'A' on the right-hand side, i.e. they are modelled as
#> exposure-affected. The reported effects are those of Zheng & van der Laan
#> (2017), identified under sequential randomization and positivity (their Lemma
#> 1). Reading them as individual-level natural effects, contrasts of each
#> subject's own counterfactual mediator, additionally requires a cross-world
#> independence assumption that is not expected to hold if such a covariate also
#> confounds the mediator-outcome relationship (Avin, Shpitser & Pearl 2005;
#> VanderWeele 2014; VanderWeele & Tchetgen Tchetgen 2017, who propose the
#> randomized interventional analogues, available here as mediation_type = "I").
#> This check reads the model formulas only and cannot verify the causal
#> structure.
```

The choice between `"I"` and `"N"` is substantive. Miles (2023)
discusses what a non-zero interventional indirect effect does and does
not establish.

### Targeted estimation

`estimator = "tmle"` replaces the simulation with the targeted minimum
loss-based estimator of Zheng & van der Laan (2017, Section 4.3), which
also handles right-censored survival outcomes. It gives Wald intervals
from the efficient influence curve, so `R` and `mc_sample` are ignored.
It accepts only lag-style recodes and no `subset`, `out_recode` or
custom models, raising an error otherwise;
[`?mediation`](https://adayim.github.io/causalMed/reference/mediation.md)
lists the requirements. Its targeted regressions are additive
main-effects working models, so the multiple robustness Zheng & van der
Laan establish for correctly specified nuisance models is not claimed
for it.

Reusing the censored data and models from above:

``` r

fit_tmle_s <- mediation(
  data           = dat_c,
  id_var         = "id",
  time_var       = "time",
  base_vars      = "V",
  exposure       = "A",
  outcome        = "Y",
  models         = models_cens,
  init_recode    = init_s,
  in_recode      = in_s,
  mediation_type = "N",
  estimator      = "tmle",
  quiet          = TRUE,
  seed           = 2026
)
#> Warning: mediation_type = "N" requested, and covariate model(s) for {L} include
#> the exposure 'A' on the right-hand side, i.e. they are modelled as
#> exposure-affected. The reported effects are those of Zheng & van der Laan
#> (2017), identified under sequential randomization and positivity (their Lemma
#> 1). Reading them as individual-level natural effects, contrasts of each
#> subject's own counterfactual mediator, additionally requires a cross-world
#> independence assumption that is not expected to hold if such a covariate also
#> confounds the mediator-outcome relationship (Avin, Shpitser & Pearl 2005;
#> VanderWeele 2014; VanderWeele & Tchetgen Tchetgen 2017, who propose the
#> randomized interventional analogues, available here as mediation_type = "I").
#> This check reads the model formulas only and cannot verify the causal
#> structure.

fit_tmle_s$estimate
#>                  Effect         RD       RR         Sd     Sd_RR   norm_lcl
#>                  <char>      <num>    <num>      <num>     <num>      <num>
#> 1:      Indirect effect  0.2525482 1.475027 0.02484973 0.0697680  0.2038436
#> 2:        Direct effect  0.2174757 1.692214 0.03729178 0.1474051  0.1443851
#> 3:         Total effect  0.4700239 2.496061 0.02964941 0.1741515  0.4119121
#> 4: Mediation Proportion 53.7309294       NA 6.07150776        NA 41.8309929
#>      norm_ucl norm_lcl_RR norm_ucl_RR
#>         <num>       <num>       <num>
#> 1:  0.3012528    1.338285    1.611770
#> 2:  0.2905662    1.403305    1.981122
#> 3:  0.5281357    2.154731    2.837392
#> 4: 65.6308660          NA          NA
```

- The TMLE changes the estimator, not the estimand, so the warning above
  still applies; [`print()`](https://rdrr.io/r/base/print.html) repeats
  it and `fit$intermediate_confounders` names the covariates.
- Where few subjects follow a regime (e.g. never exposed at every time),
  the affected targeting steps are skipped with a warning and the
  estimate relies on the fitted models; inspect these warnings. In a
  replication of the heavily censored simulation of Zheng & van der Laan
  (2017, Section 5), this implementation reproduced their value
  $`J^{1,0} \approx 0.912`$, but the intervals for the never-treated
  functionals covered about 87% rather than 95% at n = 4000, while the
  treated functionals covered at the nominal rate.

The subject-level influence-curve values are in
`fit_tmle_s$tmle_diag$eic`, for custom contrasts or diagnostics.

## References

- Avin, C., Shpitser, I., & Pearl, J. (2005). Identifiability of
  path-specific effects. *Proceedings of the 19th International Joint
  Conference on Artificial Intelligence*, 357–363.
- Keil, A. P., Edwards, J. K., Richardson, D. B., Naimi, A. I., &
  Cole, S. R. (2014). The parametric g-formula for time-to-event data:
  intuition and a worked example. *Epidemiology*, 25(6), 889–897.
- Lin, S.-H., Young, J. G., Logan, R., & VanderWeele, T. J. (2017).
  Mediation analysis for a survival outcome with time-varying exposures,
  mediators, and confounders. *Statistics in Medicine*, 36, 4153–4166.
- Miles, C. H. (2023). On the causal interpretation of randomised
  interventional indirect effects. *Journal of the Royal Statistical
  Society: Series B*, 85(4), 1154–1172.
- Robins, J. M. (1986). A new approach to causal inference in mortality
  studies with a sustained exposure period. *Mathematical Modelling*,
  7(9–12), 1393–1512.
- VanderWeele, T. J., & Tchetgen Tchetgen, E. J. (2017). Mediation
  analysis with time varying exposures and mediators. *Journal of the
  Royal Statistical Society: Series B*, 79(3), 917–938.
- Westreich, D., Cole, S. R., Young, J. G., et al. (2012). The
  parametric g-formula to estimate the effect of highly active
  antiretroviral therapy on incident AIDS or death. *Statistics in
  Medicine*, 31, 2000–2009.
- Yamamuro, S., Shinozaki, T., Iimuro, S., & Matsuyama, Y. (2021).
  Mediational g-formula for time-varying treatment and repeated-measured
  multiple mediators. *Statistical Methods in Medical Research*, 30(8),
  1782–1799.
- Zheng, W., & van der Laan, M. (2017). Longitudinal mediation analysis
  with time-varying mediators and exposures, with application to
  survival outcomes. *Journal of Causal Inference*, 5(2).
