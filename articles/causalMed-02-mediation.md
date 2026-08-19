# Causal Mediation Analysis with causalMed

## Introduction

[`mediation()`](https://adayim.github.io/causalMed/reference/mediation.md)
decomposes the total effect of a time-varying exposure into a **direct**
component (not through the mediator) and an **indirect** component
(through the mediator), for time-varying mediators and confounders —
including confounders affected by prior exposure, the setting the
mediational g-formula was developed for (Lin et al. 2017; VanderWeele &
Tchetgen Tchetgen 2017). This vignette is the full mediation story:

1.  the two estimands — interventional and natural effects — and how to
    choose,
2.  requirements and defaults,
3.  a worked single-mediator analysis on a survival outcome, reading
    every row of the output, and the `n_vw` permutation-averaging
    parameter,
4.  multiple mediators — replicating the Yamamuro et al. (2021)
    simulation, with estimates checked against documented true values,
5.  censoring,
6.  restricting models with `subset` (absorbing states),
7.  bootstrap confidence intervals, and
8.  natural effects, intermediate confounding, and the targeted (TMLE)
    estimator.

The shared vocabulary (long-format data,
[`spec_model()`](https://adayim.github.io/causalMed/reference/spec_model.md),
the
[`recodes()`](https://adayim.github.io/causalMed/reference/recodes.md)
lag hooks) is introduced in
[`vignette("causalMed-01-overview")`](https://adayim.github.io/causalMed/articles/causalMed-01-overview.md);
total-effect estimation with
[`gformula()`](https://adayim.github.io/causalMed/reference/gformula.md)
is in
[`vignette("causalMed-03-gformula")`](https://adayim.github.io/causalMed/articles/causalMed-03-gformula.md).

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

## The two estimands

For **natural** effects (`mediation_type = "N"`) the decomposition is
exact:

``` math
\underbrace{E[Y_{1,M(1)}] - E[Y_{0,M(0)}]}_{\text{Total effect}} =
  \underbrace{E[Y_{1,M(0)}] - E[Y_{0,M(0)}]}_{\text{Direct effect}} +
  \underbrace{E[Y_{1,M(1)}] - E[Y_{1,M(0)}]}_{\text{Indirect effect}}
```

where $`Y_{a, M(a^*)}`$ is the potential outcome under exposure $`a`$
with the mediator at its natural value under $`a^*`$. The mediator is
drawn from its **conditional** distribution: at each time step it is
predicted from the individual’s own covariate history, with the exposure
argument in the mediator model swapped to the other level.

For **interventional** effects (`mediation_type = "I"`, the default) the
mediator is instead set to a *stochastic* draw $`G_{a^*}`$ from its
marginal distribution under $`a^*`$ (drawn independently across
mediators). The direct and indirect effects are defined relative to
these draws and sum to the **interventional overall effect**
$`E[Y_{1,G_1}] - E[Y_{0,G_0}]`$, which generally differs from the
natural total effect $`E[Y_1] - E[Y_0]`$.
[`mediation()`](https://adayim.github.io/causalMed/reference/mediation.md)
therefore reports the total effect from a separate natural-course pass
and adds a row `TE - (Direct + Indirect)` — the mediated-interaction
residual ($`TE`$ minus the overall effect). For natural effects that
residual is exactly zero and the row is omitted.

The key differences:

|  | Interventional (`"I"`) | Natural (`"N"`) |
|----|----|----|
| Mediator drawn from | Marginal distribution (permutation) | Conditional distribution (individual history) |
| Cross-world independence required | No | Yes |
| Identified with exposure-affected confounders | Yes | **No** |
| Interpretation | Shift in mediator *population* distribution | Hypothetical individual-level swap |
| Reference | Lin et al. (2017) | Zheng & van der Laan (2017) |

The identifiability row is the practical decision point: natural effects
are **not point-identified** when a confounder of the mediator–outcome
relationship is itself affected by the exposure — an *intermediate
confounder* (Avin, Shpitser & Pearl 2005; VanderWeele & Tchetgen
Tchetgen 2017).
[`mediation()`](https://adayim.github.io/causalMed/reference/mediation.md)
warns when a covariate model carries the exposure on its right-hand
side, i.e. when a covariate is *modelled as* exposure-affected. That is
a scan of your formulas, not of the causal structure: it cannot confirm
that such a covariate also confounds the mediator–outcome relationship,
and its silence does not establish that no intermediate confounder
exists. Deciding that is yours. VanderWeele & Tchetgen Tchetgen (2017)
propose the interventional analogues for this setting, which is why
`"I"` is the default. The natural effects section below shows the
warning in action.

## Requirements and defaults

**Binary exposure.** The exposure variable must be coded as `0`
(reference/untreated) and `1` (active/treated). Both Lin et al. (2017)
and Zheng & van der Laan (2017) are defined for binary exposures;
[`mediation()`](https://adayim.github.io/causalMed/reference/mediation.md)
stops with an informative error if other values are found.

**One or more mediators.** At least one `mod_type = "mediator"` model is
required. Multiple mediators are supported for `mediation_type = "I"`
(Yamamuro et al. 2021) — list them in temporal order;
`mediation_type = "N"` is single-mediator only.

**Temporal ordering.** The list order determines the simulation
sequence. Two orderings are common:

- `A → L → M → S` — confounders at time *t* are not affected by the
  mediator at time *t* (mediator conditions on L; outcome conditions on
  A, M, L).
- `A → M → L → S` — confounders at time *t* are affected by both
  exposure and mediator (Lin et al. 2017 DAG). This is the setting where
  standard regression fails.

The function checks that exposure precedes the mediator and that the
mediator precedes the outcome, and warns if either is violated.

**Joint mediator trajectory (interventional effects).** For
`mediation_type = "I"`, a natural-course Monte Carlo cohort is simulated
under each treatment level $`a^*`$ and every individual’s full mediator
trajectory $`M(1{:}T)`$ is stored in a pool. Each intervention that
fixes a mediator to its $`a^*`$ value — **including the reference
interventions** $`\Phi_{00} = E[Y_{0,G_0}]`$ and
$`\Phi_{11} = E[Y_{1,G_1}]`$, not only the cross-world $`\Phi_{10}`$ —
permutes the relevant pool once and assigns subject $`i`$ the *entire*
trajectory of one randomly chosen pool individual (each mediator
permuted independently). This is the joint-trajectory algorithm
described by Yamamuro et al. (2021, Figure 3 step 3) and implemented by
the SAS `mGFORMULA` macro (Lin et al. 2017 eAppendix); it samples the
marginal mediator distribution targeted by Lin et al. (2017, Eq. 4) and
VanderWeele & Tchetgen Tchetgen (2017). The pool is not
survival-weighted; the full reference cohort is used at every time step
(matching both reference SAS implementations).

**Warning summary.** Warnings from model fitting (e.g., convergence,
near-separation) are held and printed as a deduplicated summary at
function exit. Repeated warnings (e.g., across 500 bootstrap replicates)
are shown once with a count.

**Outcome types.**
[`mediation()`](https://adayim.github.io/causalMed/reference/mediation.md)
handles binary end-of-follow-up outcomes (`mod_type = "outcome"`) and
survival outcomes (`mod_type = "survival"`, with optional censoring)
through the same interface — the outcome model’s `mod_type` is the only
difference. This vignette’s worked examples use a survival outcome, the
setting the mediational g-formula was developed for (Lin et al. 2017); a
compact binary-outcome example is the mediation quick start in
[`vignette("causalMed-01-overview")`](https://adayim.github.io/causalMed/articles/causalMed-01-overview.md).

------------------------------------------------------------------------

## Survival data

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
  overall effect* `Phi11 − Phi00`, not to TE; this row is the arithmetic
  difference between the two, reported so the decomposition can be read
  against the total effect. Yamamuro et al. (2021) report the
  corresponding quantity for their simulation. (For
  `mediation_type = "N"` the decomposition is exact and this row is
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

Every intervention that draws its mediators from a permuted pool does so
by randomly permuting the pool of simulated mediator trajectories. Under
`mediation_type = "I"` that is the *whole* decomposition — the
references $`\Phi_{00}`$ and $`\Phi_{11}`$ as well as the cross-world
$`\Phi_{10}`$ — but not the natural-course interventions `nat0`/`nat1`,
whose mediators come from their own fitted models. `n_vw` controls how
many independent permutations are averaged per intervention (default
`2`, matching the SAS `mGFORMULA` macro). Averaging reduces Monte Carlo
noise from the permutation step at the cost of one extra simulation pass
per intervention; `n_vw = 1` is faster and slightly noisier. It has no
effect on natural effects (`mediation_type = "N"`), which do not use
permutation.

One consequence worth knowing if you use `return_data = TRUE`: with
`n_vw > 1`, `sim_data` keeps only the last permutation of each
pool-drawing intervention while `effect_size$Est` averages all of them,
so recomputing `mean(Pred_Y)` from `sim_data` will not exactly reproduce
`Est` for those interventions. Set `n_vw = 1` if you need the two to
agree.

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
percentage points. The `True` column is a large-sample value computed
from the data-generating process, whereas the estimate comes from one
finite dataset simulated a finite number of times, so it carries both
sampling error and Monte Carlo error; the run above is a single point
estimate with neither quantified. A bootstrap (`R > 1`) is what supplies
that uncertainty in applied use.

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
censoring*, as in the g-formula treatment of right-censoring described
by Robins (1986) and Westreich et al. (2012). The censoring model’s role
is to let the hazard model be fitted on data where censoring depends on
measured covariates. Whether these estimates recover the full-data ones
depends on censoring being independent of the outcome given the modelled
covariates — an assumption about your data, not something the package
can check. Compare the two tables here to see how they came out in this
constructed example.

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
walked through in
[`vignette("causalMed-03-gformula")`](https://adayim.github.io/causalMed/articles/causalMed-03-gformula.md)
(Keil et al. 2014 describe the analysis in full, including the spline
knots).

------------------------------------------------------------------------

## Bootstrap Confidence Intervals

As elsewhere in the package, set `R > 1` for subject-level bootstrap CIs
and optionally register a parallel plan first (the general bootstrap
mechanics — percentile versus normal intervals, parallel plans — are
covered in
[`vignette("causalMed-03-gformula")`](https://adayim.github.io/causalMed/articles/causalMed-03-gformula.md)):

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
setting (see the estimands section above).
[`mediation()`](https://adayim.github.io/causalMed/reference/mediation.md)
finds the exposure on the right-hand side of a covariate model and
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
#> the exposure 'A' on the right-hand side, i.e. they are modelled as
#> exposure-affected. If such a covariate also confounds the mediator-outcome
#> relationship, natural direct and indirect effects are NOT identifiable from
#> observational data (Avin, Shpitser & Pearl 2005; VanderWeele 2014; VanderWeele
#> & Tchetgen Tchetgen 2017); for that setting VanderWeele & Tchetgen Tchetgen
#> (2017) propose the randomized interventional analogues, available here as
#> mediation_type = "I". This check reads the model formulas only and cannot
#> verify the causal structure.
```

Interventional effects (`"I"`) remain identifiable under intermediate
confounding, and VanderWeele & Tchetgen Tchetgen (2017) propose them for
exactly this setting; that is why `"I"` is the package default. The
choice between the two estimands is a substantive one — see Miles (2023)
on what a non-zero interventional indirect effect does and does not
establish.

### A targeted (TMLE) estimator for natural effects

When the natural-effects estimand *is* appropriate for your data, the
plug-in simulation is not the only estimator: `estimator = "tmle"`
implements the targeted maximum likelihood estimator from the same
reference (Zheng & van der Laan 2017, Section 4.3), including
right-censored survival outcomes. Instead of simulating forward, it runs
backward iterated regressions with targeted fluctuation steps, weighting
by inverse treatment *and* censoring probabilities. Zheng & van der Laan
(2017) establish multiple robustness for it — consistency when certain
subsets of the nuisance models are correct — and it reports Wald CIs
from the efficient influence curve, so no bootstrap is required (`R` and
`mc_sample` are ignored). Those results are stated for correctly
specified nuisance models; this implementation builds its targeted
sequential regressions as additive main-effects working models, so
transformations and interactions written into your formulas are not
carried into them (see
[`?mediation`](https://adayim.github.io/causalMed/reference/mediation.md)).

Reusing the censored dataset and model list from the censoring section
(which already includes the required exposure model):

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
#> Warning: mediation_type = "N" requested, but covariate model(s) for {L} include
#> the exposure 'A' on the right-hand side, i.e. they are modelled as
#> exposure-affected. If such a covariate also confounds the mediator-outcome
#> relationship, natural direct and indirect effects are NOT identifiable from
#> observational data (Avin, Shpitser & Pearl 2005; VanderWeele 2014; VanderWeele
#> & Tchetgen Tchetgen 2017); for that setting VanderWeele & Tchetgen Tchetgen
#> (2017) propose the randomized interventional analogues, available here as
#> mediation_type = "I". This check reads the model formulas only and cannot
#> verify the causal structure.
#> TMLE L-step t=4: only 3 observation(s) follow the intervened regime at this node (practical positivity violation). Fluctuation skipped; the estimate relies on model extrapolation here.
#>   [repeated 2 time(s)]
#> =============
#> TMLE R-step t=4: only 3 observation(s) follow the intervened regime at this node (practical positivity violation). Fluctuation skipped; the estimate relies on model extrapolation here.
#>   [repeated 2 time(s)]
#> =============
#> TMLE Z-step t=4: only 3 observation(s) follow the intervened regime at this node (practical positivity violation). Fluctuation skipped; the estimate relies on model extrapolation here.
#>   [repeated 2 time(s)]

fit_tmle_s$estimate
```

Note the **recode restriction**: the targeted engine evaluates only
lag-style recodes, so `in_recode` entries must copy a single column (as
`init_s`/`in_s` do here), lags of the exposure must copy the exposure
itself (chained lags like `lag2_A = lag_A` are rejected), `init_recode`
entries must be a constant or a column name, and `out_recode` is not
supported. Derived recodes — splines, cumulative counts, carry-forward
flags, such as the absorbing-state pattern shown earlier — and deeper
exposure history require `estimator = "gcomp"`. Violations raise an
error rather than being silently dropped, so a distorted TMLE cannot
reach you unnoticed.

Two caveats, stated plainly:

- **The identification warning still applies.** The TMLE changes the
  *estimator*, not the estimand: if an intermediate confounder makes
  natural effects unidentifiable (as in this dataset, where `L` depends
  on current `A`), no estimator can fix that. The same warning is
  emitted at run time and re-surfaced as a short caveat under the
  decomposition when you [`print()`](https://rdrr.io/r/base/print.html)
  the result (the offending confounders are also stored in
  `fit$intermediate_confounders`). It is for this setting that
  VanderWeele & Tchetgen Tchetgen (2017) propose the interventional
  estimand.
- **Positivity warnings deserve attention.** If few subjects follow an
  intervened regime (e.g., never treated at every time point), the
  affected targeting steps are skipped with a collected warning and the
  corresponding functional leans on model extrapolation — exactly as the
  plug-in does, but with the sparse support made visible instead of
  silent. Neither estimator can conjure information the data do not
  contain; the TMLE tells you when that is happening, where the plug-in
  produces a seemingly precise answer by extrapolating from its
  parametric models. When the two estimators disagree materially, the
  first things to inspect are these positivity warnings and the model
  specifications. Relatedly, the influence-curve Wald CIs can be
  somewhat anti-conservative for functionals under positivity strain: in
  our replication of the reference paper’s heavily censored simulation,
  coverage for the never-treat functionals was about 87% (rather than
  95%) even at n = 4000, while the treated functionals covered at
  nominal rates. Interpret CIs for sparsely supported regimes with
  corresponding caution.

The per-subject influence-curve values are returned in
`fit_tmle_s$tmle_diag$eic` (one column per functional) for custom
contrasts or diagnostics; `fit_tmle_s$tmle_diag$eic_mean` should be near
zero for well-supported functionals.

In simulations replicating the survival design of Zheng & van der Laan
(2017, Section 5), this implementation reproduces the published true
value of the cross-world functional and shows the expected bias
reduction from targeting, with the caveat that its GLM-based sequential
regressions are coarser than the Super Learner fits used in the paper.

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
- Avin, C., Shpitser, I., & Pearl, J. (2005). Identifiability of
  path-specific effects. *Proceedings of the 19th International Joint
  Conference on Artificial Intelligence (IJCAI)*, 357–363.
- Miles, C. H. (2023). On the causal interpretation of randomised
  interventional indirect effects. *Journal of the Royal Statistical
  Society: Series B*, 85(4), 1154–1172.
- Keil, A. P., Edwards, J. K., Richardson, D. B., Naimi, A. I., &
  Cole, S. R. (2014). The parametric g-formula for time-to-event data:
  intuition and a worked example. *Epidemiology*, 25(6), 889–897.
- Westreich, D., Cole, S. R., Young, J. G., et al. (2012). The
  parametric g-formula to estimate the effect of highly active
  antiretroviral therapy on incident AIDS or death. *Statistics in
  Medicine*, 31, 2000–2009.
