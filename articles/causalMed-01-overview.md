# Getting Started with causalMed

## Introduction

**causalMed** implements the parametric g-formula for two related goals:

1.  **Total effect estimation**
    ([`gformula()`](https://adayim.github.io/causalMed/reference/gformula.md)):
    estimate the counterfactual mean outcome had everyone in the study
    followed a given exposure strategy, using the standard parametric
    g-formula of Westreich et al. (2012) and McGrath et al. (2020).

2.  **Causal mediation analysis**
    ([`mediation()`](https://adayim.github.io/causalMed/reference/mediation.md)):
    decompose the total effect into a direct component (not through the
    mediator) and an indirect component (through the mediator), using
    either **interventional** direct/indirect effects (Lin et al. 2017)
    or **natural** direct/indirect effects (Zheng & van der Laan 2017).

Both handle time-varying exposures, mediators, and confounders —
including confounders that are themselves affected by prior exposure.
This is the setting the g-formula was introduced for (Robins 1986), and
the one in which natural direct and indirect effects are not
identifiable (Avin, Shpitser & Pearl 2005; VanderWeele & Tchetgen
Tchetgen 2017).

This vignette is the short tour: the data format, the
model-specification vocabulary that every analysis uses, and one
quick-start example of each function. The full contents live in their
own vignettes:

| Vignette | Contents |
|----|----|
| [`vignette("causalMed-02-mediation")`](https://adayim.github.io/causalMed/articles/causalMed-02-mediation.md) | Mediation: estimands, reading the decomposition, multiple mediators, censoring, natural effects, and the targeted (TMLE) estimator. |
| [`vignette("causalMed-03-gformula")`](https://adayim.github.io/causalMed/articles/causalMed-03-gformula.md) | Total effects: binary and survival outcomes, dynamic interventions, a published replication (Keil et al. 2014), bootstrap CIs, results handling, custom distributions. |
| [`vignette("causalMed-04-vs-gfoRmula")`](https://adayim.github.io/causalMed/articles/causalMed-04-vs-gfoRmula.md) | How the total-effect engine compares with the CRAN reference implementation `gfoRmula`. |

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

## Data Structure

The package requires **long-format** data: one row per subject per time
point.

``` r

data("nonsurvivaldata")
head(nonsurvivaldata, 10)
#>    id time          V         L1 L2 A          M     Y_cont Y_bin lag1_A
#> 1   1    0  0.4365731  0.7601565  0 1  0.3825810         NA    NA     NA
#> 2   1    1  0.4365731  0.1041467  1 1  1.2635708         NA    NA      1
#> 3   1    2  0.4365731  0.8956876  0 1  0.6277357         NA    NA      1
#> 4   1    3  0.4365731  1.6316564  0 1  1.2583611         NA    NA      1
#> 5   1    4  0.4365731  1.1148361  0 1  0.4602865 0.18880757     1      1
#> 6   2    0 -1.8666578 -0.3374623  0 1  0.3350378         NA    NA     NA
#> 7   2    1 -1.8666578 -0.7691214  0 1 -1.0494745         NA    NA      1
#> 8   2    2 -1.8666578 -0.7474494  0 0 -0.1559896         NA    NA      1
#> 9   2    3 -1.8666578 -1.4090944  1 0 -1.7178189         NA    NA      0
#> 10  2    4 -1.8666578 -0.3208144  1 1 -0.3123113 0.03748698     0      0
#>       lag1_L1 lag1_L2     lag1_M
#> 1          NA      NA         NA
#> 2   0.7601565       0  0.3825810
#> 3   0.1041467       1  1.2635708
#> 4   0.8956876       0  0.6277357
#> 5   1.6316564       0  1.2583611
#> 6          NA      NA         NA
#> 7  -0.3374623       0  0.3350378
#> 8  -0.7691214       0 -1.0494745
#> 9  -0.7474494       0 -0.1559896
#> 10 -1.4090944       1 -1.7178189
```

The `nonsurvivaldata` dataset contains 3 000 subjects observed at five
time points (0, 1, 2, 3, 4), 15 000 rows in total:

| Variable | Role                               |
|----------|------------------------------------|
| `id`     | Subject identifier                 |
| `time`   | Time index (0, 1, 2, 3, 4)         |
| `V`      | Time-fixed baseline covariate      |
| `A`      | Time-varying binary exposure       |
| `L1`     | Time-varying continuous confounder |
| `L2`     | Time-varying binary confounder     |
| `M`      | Time-varying continuous mediator   |
| `Y_bin`  | Binary outcome (end-of-follow-up)  |

The temporal ordering within each period is **A → L → M → Y**: the
exposure is decided first, the confounders respond to it, then the
mediator, then the outcome (this matches the data-generating process
documented in
[`?nonsurvivaldata`](https://adayim.github.io/causalMed/reference/nonsurvivaldata.md)).

------------------------------------------------------------------------

## Specifying Models

Every variable that is not time-fixed and not the identifier needs a
parametric model. Models are created with
[`spec_model()`](https://adayim.github.io/causalMed/reference/spec_model.md)
and collected into a list **in the temporal order they should be
simulated**.

``` r

spec_model(
  formula,            # Standard R formula: response ~ predictors
  var_type,           # Distribution for simulation
  mod_type,           # Role in the causal structure
  subset  = NULL,     # Optional condition restricting which rows are used
  recode  = NULL      # Optional within-loop recoding (see recodes())
)
```

### `var_type` — how to draw simulated values

The four built-in `var_type` values cover the most common cases:

| `var_type` | Distribution |
|----|----|
| `"binary"` | Bernoulli (logistic regression) |
| `"normal"` | Gaussian (linear regression), clipped to the observed range by default |
| `"categorical"` | Multinomial ([`nnet::multinom`](https://rdrr.io/pkg/nnet/man/multinom.html)) |
| `"custom"` | User-supplied fit and simulation functions |

**Simulated numeric values are clipped to the observed range of the
response by default.** `gfoRmula` applies the same rule by default,
through its `sim_trunc` argument. Pass `spec_model(truncate = FALSE)` to
draw from the untruncated fitted distribution instead (this also stops
the clipping being applied to `custom_sim` output). Which of the two is
appropriate depends on the variable being simulated.

For distributions still not covered — **truncated normal**,
**zero-inflated normal**, **absorbing states** — use
`var_type = "custom"` with the `custom_fit` and `custom_sim` arguments
to
[`spec_model()`](https://adayim.github.io/causalMed/reference/spec_model.md).
The example below implements a **zero-inflated normal** (a point mass at
zero, Gaussian otherwise — e.g. a biomarker that is exactly zero for
part of the population) as a two-part model. The object returned by
`custom_fit` can be *any* structure — here a plain list holding two fits
— as long as `custom_sim` knows how to use it:

``` r

# Part 1: is the response positive?  Part 2: its level, among the positives.
zin_fit <- function(formula, data, ...) {
  y <- model.response(model.frame(formula, data = data))
  fit_any <- glm(update(formula, I(. > 0) ~ .), family = binomial(), data = data)
  fit_pos <- lm(formula, data = data[y > 0, ])
  list(any = fit_any, pos = fit_pos)
}

zin_sim <- function(model, newdt, ...) {
  p_any  <- predict(model$any, newdata = newdt, type = "response")
  m_pos  <- predict(model$pos, newdata = newdt)
  is_pos <- rbinom(nrow(newdt), 1, p_any)
  pos    <- pmax(rnorm(nrow(newdt), m_pos, sigma(model$pos)), 0)
  is_pos * pos
}

m_zin <- spec_model(
  X ~ A + L + time,
  var_type   = "custom",
  mod_type   = "covariate",
  custom_fit = zin_fit,
  custom_sim = zin_sim,
  truncate   = FALSE   # zin_sim is authoritative over its own range
)
```

`custom_fit(formula, data, ...)` is called once during model fitting and
must return an object that `custom_sim` knows how to predict from.
`custom_sim(model, newdt, ...)` is called at each simulation step and
must return a vector of simulated values of length `nrow(newdt)`. (A
plain bounded normal needs no custom type at all — that is
`var_type = "normal"`’s default clipping behaviour, as described above)

### `mod_type` — causal role

| Value | Description |
|----|----|
| `"covariate"` | Time-varying confounder |
| `"exposure"` | Intervention variable |
| `"mediator"` | Mediator (required for [`mediation()`](https://adayim.github.io/causalMed/reference/mediation.md)) |
| `"outcome"` | Binary end-of-follow-up outcome |
| `"survival"` | Discrete-time event indicator (hazard model) |
| `"censor"` | Right-censoring indicator |

### Temporal ordering is critical

The list order determines the simulation sequence at each time step, and
it must follow your assumed data-generating mechanism. Confounders at
time *t* that are affected by the exposure at time *t* must come
**after** the exposure model and **before** the outcome model.

------------------------------------------------------------------------

## Managing Lagged Variables

Most models condition on previous-period values (lags). These are
created with
[`recodes()`](https://adayim.github.io/causalMed/reference/recodes.md),
which captures expressions to be evaluated inside the simulated dataset
at each step.

| Hook | When applied | Typical use |
|----|----|----|
| `init_recode` | Once at the first time point | Set lag variables to their baseline values |
| `in_recode` | Start of each subsequent time step | Update lags from the previous step’s values |
| `out_recode` | End of each time step | Post-simulation transforms (cumulative sums, etc.) |

``` r

init_rc <- recodes(lag1_A  = 0,   # At t=0, all lags initialised to 0
                   lag1_L1 = 0,
                   lag1_L2 = 0)

in_rc   <- recodes(lag1_A  = A,   # At each subsequent step, copy current values
                   lag1_L1 = L1,
                   lag1_L2 = L2)
```

The hooks are required even when the lag columns already exist in your
data. The Monte Carlo cohort is built from `id_var` and `base_vars`
only, so any column a model references that is neither of those — every
lag column included — has to be created by `init_recode` at the first
time point and maintained by `in_recode` thereafter. Without them the
run stops with `object 'lag1_A' not found`. Listing lag columns in
`base_vars` is not a substitute: they are not time-fixed, and
[`gformula()`](https://adayim.github.io/causalMed/reference/gformula.md)
warns that doing so distorts the sampled cohort.

------------------------------------------------------------------------

## Parallel bootstrap

The bootstrap loop uses
[`future.apply::future_lapply`](https://future.apply.futureverse.org/reference/future_lapply.html)
internally. Enable parallelism by setting a parallel plan before calling
the function:

``` r

library(future)
plan(multisession)   # use multiple R sessions in parallel

fit_par <- gformula(..., R = 500)

plan(sequential)     # restore default after use
```

On Unix/macOS, `plan(multicore)` forks the current session instead of
launching new ones; see
[`?future::plan`](https://future.futureverse.org/reference/plan.html)
for the trade-offs.

------------------------------------------------------------------------

## Quick start: a total effect with `gformula()`

Models in temporal order (**A → L1 → L2 → Y**, matching the documented
data-generating process — exposure first, current `A` in the confounder
models), then three strategies: the natural course, always treat, never
treat.

``` r

m_A  <- spec_model(A     ~ V + lag1_A + lag1_L1 + lag1_L2 + time,
                   var_type = "binary",  mod_type = "exposure")
m_L1 <- spec_model(L1    ~ V + A + lag1_L1 + time,
                   var_type = "normal",  mod_type = "covariate")
m_L2 <- spec_model(L2    ~ V + A + lag1_L2 + time,
                   var_type = "binary",  mod_type = "covariate")
m_Y  <- spec_model(Y_bin ~ V + A + L1 + L2,
                   var_type = "binary",  mod_type = "outcome")

fit_bin <- gformula(
  data        = nonsurvivaldata,
  id_var      = "id",
  time_var    = "time",
  base_vars   = "V",
  exposure    = "A",
  models      = list(m_A, m_L1, m_L2, m_Y),
  intervention = list(natural = NULL, always_treat = 1, never_treat = 0),
  ref_int     = "natural",
  init_recode = init_rc,
  in_recode   = in_rc,
  mc_sample   = 10000,
  R           = 1,       # set R > 1 for bootstrap confidence intervals
  quiet       = TRUE,
  seed        = 2025
)

fit_bin$effect_size   # mean outcome under each strategy
#>    Intervention       Est
#>          <fctr>     <num>
#> 1:      natural 0.2349247
#> 2: always_treat 0.2522388
#> 3:  never_treat 0.1066566
fit_bin$estimate      # risk difference / ratio vs the natural course
#>              Intervention  Risk_type    Estimate
#>                    <char>     <char>       <num>
#> 1: always_treat - natural Difference  0.01731407
#> 2: always_treat / natural      Ratio  1.07370051
#> 3:  never_treat - natural Difference -0.12826810
#> 4:  never_treat / natural      Ratio  0.45400334
```

That is the whole workflow: models, interventions, one call. Survival
outcomes, dynamic (rule-based) interventions, the Keil et al. (2014)
GvHD replication, bootstrap confidence intervals, extracting fitted
models and simulated data, and custom covariate distributions are all in
[`vignette("causalMed-03-gformula")`](https://adayim.github.io/causalMed/articles/causalMed-03-gformula.md).

------------------------------------------------------------------------

## Quick start: mediation with `mediation()`

Add a mediator model (`mod_type = "mediator"`) and name the outcome. The
default estimand is the **interventional** direct/indirect decomposition
(Lin et al. 2017), which stays identifiable when confounders are
affected by the exposure:

``` r

init_med <- recodes(lag1_A = 0, lag1_L1 = 0, lag1_L2 = 0, lag1_M = 0)
in_med   <- recodes(lag1_A = A, lag1_L1 = L1, lag1_L2 = L2, lag1_M = M)

models_med <- list(
  spec_model(A   ~ V + lag1_L1 + lag1_L2 + lag1_A + time,
             var_type = "binary",  mod_type = "exposure"),
  spec_model(L1  ~ V + A + lag1_L1 + time,
             var_type = "normal",  mod_type = "covariate"),
  spec_model(L2  ~ V + A + lag1_L2 + time,
             var_type = "binary",  mod_type = "covariate"),
  spec_model(M   ~ V + A + L1 + L2 + lag1_M + time,
             var_type = "normal",  mod_type = "mediator"),   # <-- mediator
  spec_model(Y_bin ~ V + A + M + L1 + L2,
             var_type = "binary",  mod_type = "outcome")
)

fit_ide <- mediation(
  data           = nonsurvivaldata,
  id_var         = "id",
  time_var       = "time",
  base_vars      = "V",
  exposure       = "A",
  outcome        = "Y_bin",
  models         = models_med,
  init_recode    = init_med,
  in_recode      = in_med,
  mediation_type = "I",     # interventional IDE/IIE (the default)
  mc_sample      = 10000,
  R              = 1,
  quiet          = TRUE,
  seed           = 2025
)

fit_ide$estimate
#>                                   Effect           RD       RR
#>                                   <char>        <num>    <num>
#> 1:                       Indirect effect  0.068268280 1.413236
#> 2:                         Direct effect  0.087924817 2.137753
#> 3:                          Total effect  0.161254808 2.678705
#> 4:              TE - (Direct + Indirect)  0.005061711       NA
#> 5:                  Mediation Proportion 45.474607398       NA
#> 6: Mediation Proportion (multiplicative) 43.707616415       NA
```

The `estimate` table is the decomposition: the direct effect, the
indirect effect (through `M`), the total effect, a mediated-interaction
residual row, and the proportion mediated. A second estimand —
**natural** direct/indirect effects (`mediation_type = "N"`), with an
optional targeted maximum likelihood estimator (`estimator = "tmle"`) —
needs stronger assumptions and is not identified when a confounder
responds to the exposure (as `L1`/`L2` do here).

What the decomposition rows mean, the choice between the two estimands,
survival outcomes, multiple mediators, censoring, and the TMLE are all
in
[`vignette("causalMed-02-mediation")`](https://adayim.github.io/causalMed/articles/causalMed-02-mediation.md).

------------------------------------------------------------------------

## Causal Assumptions

The parametric g-formula identifies these effects under the following
assumptions (Robins 1986; Westreich et al. 2012; Keil et al. 2014),
which the package cannot check and which hold or fail as a property of
your data and design, not of the code:

1.  **Consistency**: the potential outcome under the observed exposure
    history equals the observed outcome.
2.  **Positivity**: every covariate pattern that occurs under the
    intervention also occurs in the observed data (non-zero probability
    of receiving each exposure level).
3.  **Sequential exchangeability**: no unmeasured confounding of the
    exposure–outcome relationship at each time point, conditional on the
    measured past.

For **natural effects** (`mediation_type = "N"`), an additional
assumption is required:

4.  **No unmeasured exposure-induced mediator–outcome confounding**:
    there are no confounders of the mediator–outcome relationship that
    are themselves caused by prior exposure. Natural direct and indirect
    effects are not identifiable when one exists (Avin, Shpitser & Pearl
    2005; VanderWeele & Tchetgen Tchetgen 2017); the randomized
    interventional analogues targeted by `mediation_type = "I"` are
    identifiable without this assumption (VanderWeele & Tchetgen
    Tchetgen 2017).

------------------------------------------------------------------------

## References

- Robins, J. M. (1986). A new approach to causal inference in mortality
  studies with a sustained exposure period — application to control of
  the healthy worker survivor effect. *Mathematical Modelling*, 7(9–12),
  1393–1512.
- Westreich, D., Cole, S. R., Young, J. G., et al. (2012). The
  parametric g-formula to estimate the effect of highly active
  antiretroviral therapy on incident AIDS or death. *Statistics in
  Medicine*, 31, 2000–2009.
- Avin, C., Shpitser, I., & Pearl, J. (2005). Identifiability of
  path-specific effects. *Proceedings of the 19th International Joint
  Conference on Artificial Intelligence (IJCAI)*, 357–363.
- VanderWeele, T. J., & Tchetgen Tchetgen, E. J. (2017). Mediation
  analysis with time varying exposures and mediators. *Journal of the
  Royal Statistical Society: Series B*, 79(3), 917–938.
- McGrath, S., Lin, V., Zhang, Z., et al. (2020). gfoRmula: An R package
  for estimating the effects of sustained treatment strategies via the
  parametric g-formula. *Patterns*, 1, 100008.
- Lin, S.-H., Young, J. G., Logan, R., & VanderWeele, T. J. (2017).
  Mediation analysis for a survival outcome with time-varying exposures,
  mediators, and confounders. *Statistics in Medicine*, 36, 4153–4166.
- Zheng, W., & van der Laan, M. (2017). Longitudinal mediation analysis
  with time-varying mediators and exposures, with application to
  survival outcomes. *Journal of Causal Inference*, 5(2), 20160006.
- Keil, A. P., Edwards, J. K., Richardson, D. B., Naimi, A. I., &
  Cole, S. R. (2014). The parametric g-formula for time-to-event data:
  intuition and a worked example. *Epidemiology*, 25(6), 889–897.
