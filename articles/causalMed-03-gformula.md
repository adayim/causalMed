# Estimating Total Effects with the Parametric G-formula

## Introduction

[`gformula()`](https://adayim.github.io/causalMed/reference/gformula.md)
estimates the mean outcome, or risk, had everyone followed a given
exposure strategy, and so the **total effect** of one strategy against
another, with the parametric g-formula (Westreich et al. 2012; McGrath
et al. 2020): it fits a model for every time-varying variable, then
simulates each subject forward in time under each intervention.

This vignette covers binary and survival outcomes, dynamic
interventions, custom covariate distributions, a published replication
(the GvHD analysis of Keil et al. 2014), bootstrap intervals, and
working with the results. Data format,
[`spec_model()`](https://adayim.github.io/causalMed/reference/spec_model.md)
and the
[`recodes()`](https://adayim.github.io/causalMed/reference/recodes.md)
hooks are introduced in
[`vignette("causalMed-01-overview")`](https://adayim.github.io/causalMed/articles/causalMed-01-overview.md);
direct and indirect effects are in
[`vignette("causalMed-02-mediation")`](https://adayim.github.io/causalMed/articles/causalMed-02-mediation.md).

``` r

library(causalMed)
library(data.table)
```

------------------------------------------------------------------------

## Example 1: binary end-of-follow-up outcome

`nonsurvivaldata` follows 3,000 subjects over five time points. Within
each period the exposure comes first and the confounders respond to it
([`?nonsurvivaldata`](https://adayim.github.io/causalMed/reference/nonsurvivaldata.md)),
so the exposure model is listed first and the confounder models include
the current `A`. The list order should always match the assumed
data-generating process.

``` r

data("nonsurvivaldata")

# Lag bookkeeping (see the overview vignette for the recode hooks)
init_rc <- recodes(lag1_A = 0, lag1_L1 = 0, lag1_L2 = 0)
in_rc   <- recodes(lag1_A = A, lag1_L1 = L1, lag1_L2 = L2)

# ── 1. Models in temporal order: A → L1 → L2 → Y ──────────────────────────
m_A  <- spec_model(A     ~ V + lag1_A + lag1_L1 + lag1_L2 + time,
                   var_type = "binary",  mod_type = "exposure")
m_L1 <- spec_model(L1    ~ V + A + lag1_L1 + time,
                   var_type = "normal",  mod_type = "covariate")
m_L2 <- spec_model(L2    ~ V + A + lag1_L2 + time,
                   var_type = "binary",  mod_type = "covariate")
m_Y  <- spec_model(Y_bin ~ V + A + L1 + L2,
                   var_type = "binary",  mod_type = "outcome")

models_bin <- list(m_A, m_L1, m_L2, m_Y)

# ── 2. Intervention strategies ─────────────────────────────────────────────
# NULL  = natural course (draw exposure from its fitted model)
# 1 / 0 = always treat / never treat
ints <- list(natural = NULL, always_treat = 1, never_treat = 0)

# ── 3. Run the g-formula ───────────────────────────────────────────────────
fit_bin <- gformula(
  data        = nonsurvivaldata,
  id_var      = "id",
  time_var    = "time",
  base_vars   = "V",
  exposure    = "A",
  models      = models_bin,
  intervention = ints,
  ref_int     = "natural",
  init_recode = init_rc,
  in_recode   = in_rc,
  mc_sample   = 10000,
  R           = 1,        # set R > 1 for bootstrap CIs (see below)
  quiet       = TRUE,
  seed        = 2025
)
```

``` r

# Risk (mean outcome) under each strategy
fit_bin$effect_size
#>    Intervention       Est
#>          <fctr>     <num>
#> 1:      natural 0.2349247
#> 2: always_treat 0.2522388
#> 3:  never_treat 0.1066566

# Contrasts vs the reference (natural course)
fit_bin$estimate
#>              Intervention  Risk_type    Estimate
#>                    <char>     <char>       <num>
#> 1: always_treat - natural Difference  0.01731407
#> 2: always_treat / natural      Ratio  1.07370051
#> 3:  never_treat - natural Difference -0.12826810
#> 4:  never_treat / natural      Ratio  0.45400334
```

`effect_size` gives the mean outcome under each intervention; `estimate`
gives the risk difference and risk ratio against the reference.

------------------------------------------------------------------------

## Example 2: survival (time-to-event) outcome

For a survival outcome, a `mod_type = "survival"` model gives the
discrete-time hazard, which is accumulated into the cumulative incidence
$`1 - \prod_t (1 - h_t)`$: the reported quantities are risks of the
event by the end of follow-up. The data need one row per subject per
period **at risk**, with no rows after the event
([`?survivaldata`](https://adayim.github.io/causalMed/reference/survivaldata.md)).

``` r

data("survivaldata")

m_A2 <- spec_model(A ~ V + lag1_A + lag1_L + time,
                   var_type = "binary", mod_type = "exposure")
m_L  <- spec_model(L ~ V + A + lag1_L + time,
                   var_type = "normal", mod_type = "covariate")
m_Y2 <- spec_model(Y ~ V + A + L + time,
                   var_type = "binary", mod_type = "survival")  # <-- survival

models_surv <- list(m_A2, m_L, m_Y2)

fit_surv <- gformula(
  data        = survivaldata,
  id_var      = "id",
  base_vars   = "V",
  exposure    = "A",
  time_var    = "time",
  models      = models_surv,
  intervention = list(natural = NULL, never = 0, always = 1),
  ref_int     = "natural",
  init_recode = recodes(lag1_L = 0, lag1_A = 0),
  in_recode   = recodes(lag1_L = L, lag1_A = A),
  mc_sample   = 10000,
  R           = 1,
  quiet       = TRUE,
  seed        = 2025
)

fit_surv$effect_size   # Cumulative incidence by strategy
#>    Intervention       Est
#>          <fctr>     <num>
#> 1:      natural 0.6453720
#> 2:        never 0.4012772
#> 3:       always 0.7868445
fit_surv$estimate      # Risk contrasts
#>        Intervention  Risk_type   Estimate
#>              <char>     <char>      <num>
#> 1:  never - natural Difference -0.2440949
#> 2:  never / natural      Ratio  0.6217765
#> 3: always - natural Difference  0.1414724
#> 4: always / natural      Ratio  1.2192106
```

------------------------------------------------------------------------

## Dynamic (threshold) interventions

A **dynamic intervention** sets the exposure by a rule applied to each
subject’s simulated values, written with
[`dyn_int()`](https://adayim.github.io/causalMed/reference/dyn_int.md)
and evaluated in the simulated data at every time step.

The rule runs at the exposure model’s position in the list. Variables
listed *earlier* (and the exposure’s own natural-course draw) hold
current-period values; variables listed *later*, here `L1`, still hold
the previous period’s values. A rule for an exposure decided at the
start of each period therefore uses the previous period’s covariates:

``` r

fit_dyn <- gformula(
  data        = nonsurvivaldata,
  id_var      = "id",
  time_var    = "time",
  base_vars   = "V",
  exposure    = "A",
  models      = models_bin,
  intervention = list(
    natural = NULL,
    treat_if_prev_L1_pos = dyn_int(as.numeric(lag1_L1 > 0))
  ),
  ref_int     = "natural",
  init_recode = init_rc,
  in_recode   = in_rc,
  mc_sample   = 10000,
  R           = 1,
  quiet       = TRUE,
  seed        = 2025
)

fit_dyn$effect_size
#>            Intervention       Est
#>                  <fctr>     <num>
#> 1:              natural 0.2349247
#> 2: treat_if_prev_L1_pos 0.2223531
fit_dyn$estimate
#>                      Intervention  Risk_type    Estimate
#>                            <char>     <char>       <num>
#> 1: treat_if_prev_L1_pos - natural Difference -0.01257164
#> 2: treat_if_prev_L1_pos / natural      Ratio  0.94648651
```

Any expression over the columns in scope is allowed, e.g.
`dyn_int(as.numeric(A > 0 & lag1_L1 > median(lag1_L1)))`, where `A` is
the natural-course draw that the rule overrides.

------------------------------------------------------------------------

## Custom covariate distributions

`var_type = "custom"` covers distributions without a built-in type, such
as a zero-inflated normal (a point mass at zero, Gaussian otherwise),
here as a two-part model. `custom_fit` may return any object that
`custom_sim` knows how to use, here a list of two fits:

``` r

# Part 1: is the value positive?  Part 2: its level, among the positives.
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
  truncate   = FALSE   # zin_sim controls its own range
)
```

`custom_fit(formula, data, ...)` is called once, when the models are
fitted; `custom_sim(model, newdata)` is called at every time step and
returns one value per row. `custom_fit` must be defined at the top level
(or be a namespace-qualified function) so that it can be found when the
models are fitted. Choosing and checking the distribution is the
analyst’s responsibility.

------------------------------------------------------------------------

## A published example: preventing GvHD (Keil et al. 2014)

The package ships `gvhd`, the person-day bone marrow transplant data of
the parametric g-formula illustration of Keil et al. (2014). The code
below reproduces that analysis: the risk of death by day 1825 had
graft-versus-host disease (GvHD) **never** occurred, against the natural
course. It combines several features:

- **five models**: relapse → platelet recovery → GvHD (exposure) →
  censoring → death (hazard). Here the covariates are measured before
  the day’s exposure, so they come first, as in the paper;
- **absorbing states**: each state is modelled only among those not yet
  in it (`subset =`), and `out_recode` keeps it at 1 afterwards;
- a **censoring model**;
- **restricted cubic splines** of age and day, and day counters, built
  with the three recode hooks;
- daily data: 137 subjects over 1,825 days.

The models follow Appendix 2 of the paper
([`?gvhd`](https://adayim.github.io/causalMed/reference/gvhd.md)). First
the time-fixed transforms: restricted cubic splines of age (knots 17,
25.4, 30, 41.4) and of day (knots 83.6, 401.4, 947, 1862.2), and their
squares:

``` r

data("gvhd")

gvhd <- within(gvhd, {
  agesq    <- age^2
  agecurs1 <- (age > 17.0) * (age - 17.0)^3 -
              ((age > 30.0) * (age - 30.0)^3) * (41.4 - 17.0) / (41.4 - 30.0)
  agecurs2 <- (age > 25.4) * (age - 25.4)^3 -
              ((age > 41.4) * (age - 41.4)^3) * (41.4 - 25.4) / (41.4 - 30.0)
  daysq    <- day^2
  daycurs1 <- (day > 83.6)   * ((day - 83.6)   / 83.6)^3 +
              (day > 1862.2) * ((day - 1862.2) / 83.6)^3 * (947.0 - 83.6) -
              (day > 947.0)  * ((day - 947.0)  / 83.6)^3 * (1862.2 - 83.6) / (1862.2 - 947.0)
  daycurs2 <- (day > 401.4)  * ((day - 401.4)  / 83.6)^3 +
              (day > 1862.2) * ((day - 1862.2) / 83.6)^3 * (947.0 - 401.4) -
              (day > 947.0)  * ((day - 947.0)  / 83.6)^3 * (1862.2 - 401.4) / (1862.2 - 947.0)
})
```

The five models. Each state’s model applies only among those not yet in
it (`subset = ...m1 == 0`), and the death hazard interacts the day
spline with `gvhd`:

``` r

models_gvhd <- list(
  spec_model(relapse ~ all + cmv + male + age + gvhdm1 + daysgvhd + platnormm1 +
               daysnoplatnorm + agecurs1 + agecurs2 + day + daysq + wait,
             var_type = "binary", mod_type = "covariate", subset = relapsem1 == 0),
  spec_model(platnorm ~ all + cmv + male + age + agecurs1 + agecurs2 + gvhdm1 +
               daysgvhd + daysnorelapse + wait,
             var_type = "binary", mod_type = "covariate", subset = platnormm1 == 0),
  spec_model(gvhd ~ all + cmv + male + age + platnormm1 + daysnoplatnorm +
               relapsem1 + daysnorelapse + agecurs1 + agecurs2 + day + daysq + wait,
             var_type = "binary", mod_type = "exposure", subset = gvhdm1 == 0),
  spec_model(censlost ~ all + cmv + male + age + daysgvhd + daysnoplatnorm +
               daysnorelapse + agesq + day + daycurs1 + daycurs2 + wait,
             var_type = "binary", mod_type = "censor"),
  spec_model(d ~ all + cmv + male + age + gvhd + platnorm + daysnoplatnorm +
               relapse + daysnorelapse + agesq + wait +
               day * gvhd + daycurs1 * gvhd + daycurs2 * gvhd,
             var_type = "binary", mod_type = "survival")
)
```

The recode hooks: `init_recode` sets day 1 (states and counters at 0,
functions of day computed), `in_recode` updates the functions of day and
the one-day lags at the start of each later day, and `out_recode`
advances the counters and carries the absorbing states forward at the
end of each day:

``` r

init_recode <- recodes(
  daysq    = day^2,
  daycurs1 = (day > 83.6)   * ((day - 83.6)   / 83.6)^3 +
             (day > 1862.2) * ((day - 1862.2) / 83.6)^3 * (947.0 - 83.6) -
             (day > 947.0)  * ((day - 947.0)  / 83.6)^3 * (1862.2 - 83.6) / (1862.2 - 947.0),
  daycurs2 = (day > 401.4)  * ((day - 401.4)  / 83.6)^3 +
             (day > 1862.2) * ((day - 1862.2) / 83.6)^3 * (947.0 - 401.4) -
             (day > 947.0)  * ((day - 947.0)  / 83.6)^3 * (1862.2 - 401.4) / (1862.2 - 947.0),
  relapse = 0, gvhd = 0, platnorm = 0, gvhdm1 = 0, relapsem1 = 0, platnormm1 = 0,
  daysnorelapse = 0, daysnoplatnorm = 0, daysnogvhd = 0,
  daysrelapse = 0, daysplatnorm = 0, daysgvhd = 0)

in_recode <- recodes(
  daysq    = day^2,
  daycurs1 = (day > 83.6)   * ((day - 83.6)   / 83.6)^3 +
             (day > 1862.2) * ((day - 1862.2) / 83.6)^3 * (947.0 - 83.6) -
             (day > 947.0)  * ((day - 947.0)  / 83.6)^3 * (1862.2 - 83.6) / (1862.2 - 947.0),
  daycurs2 = (day > 401.4)  * ((day - 401.4)  / 83.6)^3 +
             (day > 1862.2) * ((day - 1862.2) / 83.6)^3 * (947.0 - 401.4) -
             (day > 947.0)  * ((day - 947.0)  / 83.6)^3 * (1862.2 - 401.4) / (1862.2 - 947.0),
  platnormm1 = platnorm, relapsem1 = relapse, gvhdm1 = gvhd)

out_recode <- recodes(
  daysnorelapse  = ifelse(relapse == 0,  daysnorelapse + 1,  daysnorelapse),
  daysrelapse    = ifelse(relapse == 1,  daysrelapse + 1,    daysrelapse),
  daysnoplatnorm = ifelse(platnorm == 0, daysnoplatnorm + 1, daysnoplatnorm),
  daysplatnorm   = ifelse(platnorm == 1, daysplatnorm + 1,   daysplatnorm),
  daysnogvhd     = ifelse(gvhd == 0,     daysnogvhd + 1,     daysnogvhd),
  daysgvhd       = ifelse(gvhd == 1,     daysgvhd + 1,       daysgvhd),
  # absorbing carry-forward: once a state was 1 yesterday, keep it at 1
  platnorm = ifelse(platnormm1 == 1, 1, platnorm),
  relapse  = ifelse(relapsem1 == 1,  1, relapse),
  gvhd     = ifelse(gvhdm1 == 1,     1, gvhd))
```

``` r

fit_gvhd <- gformula(gvhd,
  id_var    = "id", time_var = "day", exposure = "gvhd",
  base_vars = c("age", "agesq", "agecurs1", "agecurs2", "male", "cmv", "all", "wait"),
  models    = models_gvhd,
  intervention = list(never = 0),     # a natural-course reference is added automatically
  init_recode = init_recode, in_recode = in_recode, out_recode = out_recode,
  mc_sample = 20000, R = 1, quiet = TRUE, seed = 20260703)

fit_gvhd$effect_size
#>    Intervention       Est
#>          <fctr>     <num>
#> 1:      natural 0.6093279
#> 2:        never 0.5881213
fit_gvhd$estimate
#>       Intervention  Risk_type    Estimate
#>             <char>     <char>       <num>
#> 1: never - natural Difference -0.02120662
#> 2: never / natural      Ratio  0.96519671
```

These are the simulated 5-year risks of death under the natural course
and under “never GvHD”, and their contrast, following the specification
of Keil et al. (2014); see that paper for the interpretation and the
assumptions. With 1,825 daily time steps the simulation takes a few
minutes, and bootstrap intervals (`R > 1`) take much longer.

------------------------------------------------------------------------

## Bootstrap Confidence Intervals

Set `R > 1` for percentile and normal-approximation confidence
intervals. The bootstrap resamples whole subjects, keeping each
subject’s time points together.

``` r

fit_boot <- gformula(
  data        = nonsurvivaldata,
  id_var      = "id",
  time_var    = "time",
  base_vars   = "V",
  exposure    = "A",
  models      = models_bin,
  intervention = list(natural = NULL, always = 1),
  ref_int     = "natural",
  init_recode = init_rc,
  in_recode   = in_rc,
  mc_sample   = 10000,
  R           = 200,       # 200 bootstrap replicates
  quiet       = TRUE,
  seed        = 2025
)

# effect_size now includes Sd, perct_lcl/ucl, norm_lcl/ucl
fit_boot$effect_size
#>    Intervention       Est          Sd perct_lcl perct_ucl  norm_lcl  norm_ucl
#>          <fctr>     <num>       <num>     <num>     <num>     <num>     <num>
#> 1:      natural 0.2349247 0.007879589 0.2182889 0.2468377 0.2194810 0.2503684
#> 2:       always 0.2522388 0.008566108 0.2358453 0.2669507 0.2354495 0.2690280
fit_boot$estimate
#>        Intervention  Risk_type   Estimate          Sd  perct_lcl perct_ucl
#>              <char>     <char>      <num>       <num>      <num>     <num>
#> 1: always - natural Difference 0.01731407 0.002329456 0.01294381 0.0224198
#> 2: always / natural      Ratio 1.07370051 0.009868665 1.05635827 1.0949960
#>      norm_lcl   norm_ucl
#>         <num>      <num>
#> 1: 0.01274842 0.02187972
#> 2: 1.05435829 1.09304274

# the individual per-replicate draws are retained in boot_estimates
# ($interventions and $contrasts), for custom intervals or diagnostics
head(fit_boot$boot_estimates$interventions)
#>    replicate Intervention       Est
#>        <int>       <fctr>     <num>
#> 1:         1      natural 0.2373405
#> 2:         1       always 0.2543564
#> 3:         2      natural 0.2299337
#> 4:         2       always 0.2511462
#> 5:         3      natural 0.2310667
#> 6:         3       always 0.2498870
```

To run the replicates in parallel, set a plan first, e.g.
`future::plan(future::multisession)`.

------------------------------------------------------------------------

## Working with Results

### Extracting fitted models

`return_fitted = TRUE` returns the full fitted model objects:

``` r

fit_full <- gformula(
  data        = nonsurvivaldata,
  id_var      = "id",
  time_var    = "time",
  base_vars   = "V",
  exposure    = "A",
  models      = models_bin,
  intervention = list(natural = NULL, always = 1),
  init_recode = init_rc,
  in_recode   = in_rc,
  mc_sample   = 5000,
  R           = 1,
  return_fitted = TRUE,
  quiet       = TRUE,
  seed        = 2025
)

# Names correspond to the response variable of each model
names(fit_full$fitted_models)
#> [1] "A"     "L1"    "L2"    "Y_bin"

# Access a specific model's coefficients
coef(fit_full$fitted_models$A)
#> (Intercept)           V      lag1_A     lag1_L1     lag1_L2        time 
#>  1.29379895  0.53935295  0.17102958  0.20813430  0.25822156  0.01266473
```

### Retrieving the simulated data

`return_data = TRUE` returns the simulated data (it can be large) as an
**end-of-follow-up snapshot**: one row per simulated subject per
intervention, each variable at its last simulated time step, with the
accumulated `Pred_Y`. Earlier time steps are not kept.

``` r

fit_data <- gformula(..., return_data = TRUE)

# One row per MC subject per intervention, at the last time point
head(fit_data$sim_data)
```

------------------------------------------------------------------------

## References

- Westreich, D., Cole, S. R., Young, J. G., et al. (2012). The
  parametric g-formula to estimate the effect of highly active
  antiretroviral therapy on incident AIDS or death. *Statistics in
  Medicine*, 31, 2000–2009.
- Keil, A. P., Edwards, J. K., Richardson, D. B., Naimi, A. I., &
  Cole, S. R. (2014). The parametric g-formula for time-to-event data:
  intuition and a worked example. *Epidemiology*, 25(6), 889–897.
- McGrath, S., Lin, V., Zhang, Z., et al. (2020). gfoRmula: An R package
  for estimating the effects of sustained treatment strategies via the
  parametric g-formula. *Patterns*, 1, 100008.
- Robins, J. M. (1986). A new approach to causal inference in mortality
  studies with a sustained exposure period. *Mathematical Modelling*,
  7(9–12), 1393–1512.
