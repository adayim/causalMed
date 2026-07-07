# causalMed and gfoRmula: A Comparative Overview

## Overview

`gfoRmula` (McGrath et al. 2020) is the established CRAN reference
implementation of the parametric g-formula in R. `causalMed` builds
directly on the same statistical foundation and treats `gfoRmula` as its
benchmark for total-effect estimation.

This vignette describes:

1.  The shared statistical foundation of both packages
2.  Their common lineage in the GFORMULA-SAS reference macro
3.  What `gfoRmula` offers that `causalMed` currently does not
4.  The mediation extension that motivated `causalMed`
5.  Side-by-side total-effect estimates confirming numerical equivalence
6.  Guidance on when to use each package

``` r

library(causalMed)
library(data.table)
```

------------------------------------------------------------------------

## Shared Statistical Foundation

Both packages implement the same parametric g-formula algorithm: Monte
Carlo forward simulation under user-specified interventions on a
longitudinal dataset in long format (one row per subject per time
point). At each time step, covariate and outcome models fitted on the
observed data are applied in temporal order, the intervention is imposed
on the exposure variable, and the predicted outcomes are accumulated
across simulated individuals. Confidence intervals are obtained by
non-parametric bootstrap resampling at the individual level.

------------------------------------------------------------------------

## Relationship to GFORMULA-SAS

The original reference implementation of the parametric g-formula is the
**GFORMULA SAS macro** (Logan, Young, Taubman, Hernán and colleagues;
available from the Harvard Causal Inference group). `gfoRmula` is the R
port of that macro, and `causalMed`’s total-effect engine
([`gformula()`](https://adayim.github.io/causalMed/reference/gformula.md))
targets the same algorithm. All three share the identical core: fit
parametric models on the observed long data, Monte Carlo–resample
baseline records, forward-simulate the covariate → exposure → outcome
sequence one time step at a time, compute survival risk as
$`1 - \prod_t (1 - h_t)`$ with censoring disabled under intervention,
and bootstrap over subjects for confidence intervals. The `causalMed`
regression tests cross-check total-effect estimates against `gfoRmula`
(the R twin of the SAS macro) to Monte Carlo noise.

Two clarifications are worth keeping in mind:

- **GFORMULA-SAS estimates total effects only.** It contains no mediator
  decomposition. The direct/indirect machinery in `causalMed` (see the
  mediation section below) is an *extension* beyond what GFORMULA-SAS or
  `gfoRmula` provide. The cross-world mediator-permutation step
  `causalMed` uses for interventional effects follows a *separate* SAS
  macro, `mGFORMULA` (Lin et al. 2017), not the total-effect GFORMULA
  macro.
- **The SAS macro has a much larger declarative vocabulary.** It exposes
  up to eight covariate distribution types, keyword-driven functional
  forms of covariate history (lags, cumulative averages, restricted
  cubic splines via `ptype`), several packaged intervention types
  (threshold, increment, sampling from the observed distribution),
  competing-event handling, and additional end-of-follow-up outcome
  models. `causalMed` covers the same modelling *power* through a
  smaller surface — model formulas plus
  `recode`/[`dyn_int()`](https://adayim.github.io/causalMed/reference/dyn_int.md)
  expressions and `custom_fit`/`custom_sim` — but several of these
  conveniences must be coded by hand rather than requested by keyword.
  The specific gaps are itemised below.

------------------------------------------------------------------------

## What gfoRmula Offers That causalMed Does Not

`gfoRmula` is the more mature and feature-complete package for
total-effect g-formula estimation. Users with purely total-effect
analyses should consider `gfoRmula` as their primary tool. Key
capabilities in `gfoRmula` not currently available in `causalMed`
include:

### Richer covariate distribution types

`gfoRmula` supports several covariate types beyond binary and normal
that `causalMed` does not yet implement:

| Covariate type                                       | gfoRmula | causalMed |
|------------------------------------------------------|:--------:|:---------:|
| `"binary"`                                           |    ✓     |     ✓     |
| `"normal"`                                           |    ✓     |     ✓     |
| `"categorical"`                                      |    ✓     |     ✓     |
| `"bounded normal"` (truncated to observed range)     |    ✓     |     —     |
| `"zero-inflated normal"` (point mass at 0)           |    ✓     |     —     |
| `"truncated normal"` (left-truncated)                |    ✓     |     —     |
| `"absorbing"` (once 1, always 1)                     |    ✓     |     —     |
| `"categorical time"` (time as categorical predictor) |    ✓     |     —     |

For datasets containing variables that require these distributions
(e.g. a biomarker that can be zero or a treatment that is absorbing once
initiated), `gfoRmula` is the appropriate choice. However, one can use
`custom_fit` and `custom_sim` to implement custom distributions and
simulations in `spec_model` function. This gives users the flexibility
to model complex scenarios, but it should be done with caution and
careful validation, as it requires more manual coding and is not yet
supported by built-in diagnostics.

### Competing events

`gfoRmula` provides explicit support for competing events via the
`compevent_name` argument, modelling and eliminating the competing risk
in the simulation. `causalMed` handles censoring but does not have a
dedicated competing-event interface.

### Full-function custom interventions

`gfoRmula`’s custom intervention interface accepts arbitrary R functions
with access to the full simulated dataset, the time index, and all
parameter values at each step. This gives more control for complex
multi-variable interventions. `causalMed` offers
[`dyn_int()`](https://adayim.github.io/causalMed/reference/dyn_int.md),
which captures an R expression evaluated inside the simulated dataset at
each step. The current time step is in scope — the `time_var` column is
refreshed before every step — so a rule may reference it directly
(e.g. `dyn_int(as.numeric(A > 0 & time >= 2))`), and history-dependent
rules can be built by maintaining lagged columns through the recode
mechanism (`in_recode` / `out_recode`). The practical difference is that
`causalMed` does not pass a single pooled data-frame object to a user
function the way `gfoRmula`’s full function interface does; multi-step
logic is expressed instead through recodes plus the captured expression.

### CRAN stability and community support

`gfoRmula` is an established CRAN package with comprehensive
documentation, a dedicated publication (McGrath et al. 2020), and
ongoing development by the CausalAB group at Harvard. `causalMed` is a
development-stage package and should be used with appropriate caution.

------------------------------------------------------------------------

## What causalMed Adds: Causal Mediation Analysis

The sole motivation for `causalMed` is to extend the standard parametric
g-formula with the **survival mediational g-formula** (Lin et al. 2017;
Zheng & van der Laan 2017) for decomposing total effects into direct and
indirect components. `gfoRmula` does not support mediation analysis.

`causalMed` offers two mediation estimands, selected via
`mediation_type`:

- **`"I"` — Interventional IDE/IIE** (Lin et al. 2017): the mediator
  distribution is marginalised over confounders by randomly permuting
  mediator values simulated under the reference exposure. Does not
  require the cross-world independence assumption.

- **`"N"` — Natural NDE/NIE** (Zheng & van der Laan 2017): the mediator
  model is evaluated at the alternative exposure level while keeping the
  individual’s own covariate history. Requires stronger sequential
  no-unmeasured-confounding assumptions.

``` r

# Model list must include a mediator (mod_type = "mediator")
models_med <- list(
  spec_model(L ~ V + A_lag1 + L_lag1 + time,
             var_type = "normal", mod_type = "covariate"),
  spec_model(A ~ V + A_lag1 + L + time,
             var_type = "binary", mod_type = "exposure"),
  spec_model(M ~ V + A + L + M_lag1 + time,
             var_type = "normal", mod_type = "mediator"),
  spec_model(Y ~ V + A + M + L,
             var_type = "binary", mod_type = "outcome")
)

fit_med <- mediation(
  data           = dat_med,
  id_var         = "id",
  time_var       = "time",
  base_vars      = "V",
  exposure       = "A",
  outcome        = "Y",
  models         = models_med,
  mediation_type = "I",     # interventional IDE/IIE (Lin et al. 2017)
  init_recode    = recodes(A_lag1 = 0, L_lag1 = 0, M_lag1 = 0),
  in_recode      = recodes(A_lag1 = A, L_lag1 = L, M_lag1 = M),
  mc_sample      = 10000,
  R              = 200,
  seed           = 20250915
)

fit_med$estimate
```

The `estimate` table returns the indirect effect, direct effect, total
effect, and proportion mediated on both the risk-difference and
risk-ratio scales. See
[`vignette("causalMed-overview")`](https://adayim.github.io/causalMed/articles/causalMed-overview.md)
for a full worked example and
[`vignette("mediation-survival")`](https://adayim.github.io/causalMed/articles/mediation-survival.md)
for mediation with survival outcomes, multiple mediators, and censoring.

------------------------------------------------------------------------

## API Differences

Beyond features, the two packages differ in how models and history
functions are specified:

| Aspect | causalMed | gfoRmula |
|----|----|----|
| Model specification | [`spec_model()`](https://adayim.github.io/causalMed/reference/spec_model.md) per variable: formula + `var_type` + `mod_type` in a list | Separate `covparams`, `ymodel`, `covtypes`, `covnames`, `outcome_type` arguments |
| Lag / history | `init_recode` / `in_recode` via [`recodes()`](https://adayim.github.io/causalMed/reference/recodes.md) (arbitrary expressions) | Built-in `lagged`, `cumavg` functions via `histories` / `histvars` |
| End-of-step transforms | `out_recode` hook | Not available |
| Interventions | Named list, any number | Up to 3 (`intervention1.X`, `intervention2.X`, `intervention3.X`) |
| Risk contrasts | Computed automatically (`ref_int`) | Manual post-processing of `$result` |
| Output | S3 object with `print`/`summary` | List with `$result` data.table |

These are design choices reflecting different trade-offs, not
deficiencies in either package. `gfoRmula`’s explicit argument structure
makes each model’s role transparent; `causalMed`’s list-based approach
may be more concise when many models are specified.

------------------------------------------------------------------------

## Numerical Validation: Side-by-Side Examples

For total-effect analyses the two packages should produce equivalent
results. We illustrate this with both a binary end-of-follow-up outcome
and a survival outcome.

### Data preparation

``` r

data("nonsurvivaldata", package = "causalMed")
dat <- as.data.table(nonsurvivaldata)
dat[, time := as.integer(time)]
setorder(dat, id, time)

# Lag recodes (applied inside the Monte Carlo loop)
init_rc <- recodes(lag1_A = 0, lag1_L1 = 0, lag1_L2 = 0)
in_rc   <- recodes(lag1_A = A, lag1_L1 = L1, lag1_L2 = L2)
```

### Example 1: Binary end-of-follow-up outcome

**causalMed**

``` r

models_cm <- list(
  spec_model(L1 ~ lag1_A + lag1_L1 + V + time,
             var_type = "normal", mod_type = "covariate"),
  spec_model(L2 ~ lag1_A + lag1_L2 + V + time,
             var_type = "binary", mod_type = "covariate"),
  spec_model(A  ~ lag1_A + L1 + L2 + V + time,
             var_type = "binary", mod_type = "exposure"),
  spec_model(Y_bin ~ A + L1 + L2,
             var_type = "binary", mod_type = "outcome")
)

fit_cm <- causalMed::gformula(
  data         = dat,
  id_var       = "id",
  time_var     = "time",
  base_vars    = "V",
  exposure     = "A",
  models       = models_cm,
  intervention = list(natural = NULL, always = 1, never = 0),
  ref_int      = "natural",
  init_recode  = init_rc,
  in_recode    = in_rc,
  mc_sample    = 10000,
  R            = 1,
  quiet        = TRUE,
  seed         = 20250915
)

fit_cm$effect_size
#>    Intervention        Est
#>          <fctr>      <num>
#> 1:      natural 0.14062380
#> 2:       always 0.15087156
#> 3:        never 0.08641565
```

**gfoRmula**

``` r

suppressPackageStartupMessages(requireNamespace("gfoRmula"))

fit_gf <- gfoRmula::gformula(
  obs_data     = copy(dat),
  id           = "id",
  time_name    = "time",
  time_points  = length(unique(dat$time)),
  covnames     = c("L1", "L2", "A"),
  covtypes     = c("normal", "binary", "binary"),
  basecovs     = "V",
  covparams    = list(covmodels = c(
    L1 ~ lag1_A + lag1_L1 + V + time,
    L2 ~ lag1_A + lag1_L2 + V + time,
    A  ~ lag1_A + L1 + L2 + V + time
  )),
  ymodel       = Y_bin ~ A + L1 + L2,
  outcome_name = "Y_bin",
  outcome_type = "binary_eof",
  histories    = c(gfoRmula::lagged),
  histvars     = list(c("A", "L1", "L2")),
  intervention1.A = list(gfoRmula::static, rep(1, length(unique(dat$time)))),
  intervention2.A = list(gfoRmula::static, rep(0, length(unique(dat$time)))),
  int_descript = c("Always treat", "Never treat"),
  nsimul       = 10000,
  seed         = 20250915
)

fit_gf$result[, c("k", "g-form mean")]
#>        k g-form mean
#>    <num>       <num>
#> 1:     4  0.14107172
#> 2:     4  0.15071907
#> 3:     4  0.08655451
```

**Side-by-side at final time point**

``` r

gf_res <- fit_gf$result[k == max(fit_gf$result$k),
                         c("Interv.", "g-form mean")]
gf_res[, Intervention := c("natural", "always", "never")]
setnames(gf_res, "g-form mean", "gfoRmula")
gf_res[, gfoRmula := round(gfoRmula, 4)]

cm_res <- fit_cm$effect_size[, .(Intervention, causalMed = round(Est, 4))]
merge(cm_res, gf_res, by = "Intervention")
#> Key: <Intervention>
#>    Intervention causalMed Interv. gfoRmula
#>          <char>     <num>   <num>    <num>
#> 1:       always    0.1509       1   0.1507
#> 2:      natural    0.1406       0   0.1411
#> 3:        never    0.0864       2   0.0866
```

------------------------------------------------------------------------

### Example 2: Survival (time-to-event) outcome

**causalMed**

``` r

data("survivaldata", package = "causalMed")
dat_s <- as.data.table(survivaldata)
setorder(dat_s, id, time)

models_surv <- list(
  spec_model(L ~ V + lag1_L + time,
             var_type = "normal", mod_type = "covariate"),
  spec_model(A ~ V + lag1_A + lag1_L + L + time,
             var_type = "binary", mod_type = "exposure"),
  spec_model(Y ~ lag1_A + A + L + lag1_L + time,
             var_type = "binary", mod_type = "survival")
)

fit_cm_s <- causalMed::gformula(
  data         = dat_s,
  id_var       = "id",
  time_var     = "time",
  base_vars    = "V",
  exposure     = "A",
  models       = models_surv,
  intervention = list(natural = NULL, never = 0, always = 1),
  ref_int      = "natural",
  init_recode  = recodes(lag1_L = 0, lag1_A = 0),
  in_recode    = recodes(lag1_L = L, lag1_A = A),
  mc_sample    = 10000,
  R            = 1,
  quiet        = TRUE,
  seed         = 20250915
)

fit_cm_s$effect_size
#>    Intervention       Est
#>          <fctr>     <num>
#> 1:      natural 0.9340304
#> 2:        never 0.8538016
#> 3:       always 0.9702749
```

**gfoRmula**

``` r

T_max <- length(unique(dat_s$time))

fit_gf_s <- gfoRmula::gformula(
  obs_data     = copy(dat_s),
  id           = "id",
  time_name    = "time",
  time_points  = T_max,
  covnames     = c("A", "L"),
  covtypes     = c("binary", "normal"),
  basecovs     = "V",
  covparams    = list(covmodels = c(
    A ~ V + lag1_A + lag1_L + L + time,
    L ~ V + lag1_L + time
  )),
  ymodel       = Y ~ lag1_A + A + L + lag1_L + time,
  outcome_name = "Y",
  outcome_type = "survival",
  histories    = c(gfoRmula::lagged),
  histvars     = list(c("A", "L")),
  intervention1.A = list(gfoRmula::static, rep(0, T_max)),
  intervention2.A = list(gfoRmula::static, rep(1, T_max)),
  int_descript = c("Never treat", "Always treat"),
  nsimul       = 10000,
  seed         = 20250915,
  nsamples     = 0
)

fit_gf_s$result[k == max(fit_gf_s$result$k), c("Interv.", "g-form risk")]
#>    Interv. g-form risk
#>      <num>       <num>
#> 1:       0   0.9090761
#> 2:       1   0.8541680
#> 3:       2   0.9703732
```

**Side-by-side at final time point**

``` r

gf_s_res <- fit_gf_s$result[k == max(fit_gf_s$result$k),
                              c("Interv.", "g-form risk")]
gf_s_res[, Intervention := c("natural", "never", "always")]
setnames(gf_s_res, "g-form risk", "gfoRmula")
gf_s_res[, gfoRmula := round(gfoRmula, 4)]

cm_s_res <- fit_cm_s$effect_size[, .(Intervention,
                                      causalMed = round(Est, 4))]
merge(cm_s_res, gf_s_res, by = "Intervention")
#> Key: <Intervention>
#>    Intervention causalMed Interv. gfoRmula
#>          <char>     <num>   <num>    <num>
#> 1:       always    0.9703       2   0.9704
#> 2:      natural    0.9340       0   0.9091
#> 3:        never    0.8538       1   0.8542
```

------------------------------------------------------------------------

### Example 3: Dynamic intervention

Both packages support dynamic (rule-based) interventions. Here we
estimate the risk under “treat only if L1 \> 0”. The syntax differs
between packages but the estimand is the same.

**causalMed** —
[`dyn_int()`](https://adayim.github.io/causalMed/reference/dyn_int.md)
captures the rule as an unevaluated expression evaluated inside the
simulated dataset at each time step:

``` r

fit_dyn <- causalMed::gformula(
  data         = dat,
  id_var       = "id",
  time_var     = "time",
  base_vars    = "V",
  exposure     = "A",
  models       = models_cm,
  intervention = list(natural = NULL, treat_if_L1_pos = dyn_int(as.numeric(L1 > 0))),
  ref_int      = "natural",
  init_recode  = init_rc,
  in_recode    = in_rc,
  mc_sample    = 5000,
  R            = 1,
  quiet        = TRUE,
  seed         = 20250915
)

fit_dyn$effect_size
#>       Intervention       Est
#>             <fctr>     <num>
#> 1:         natural 0.1414215
#> 2: treat_if_L1_pos 0.1324326
fit_dyn$estimate
#>                 Intervention  Risk_type    Estimate
#>                       <char>     <char>       <num>
#> 1: treat_if_L1_pos - natural Difference -0.00898884
#> 2: treat_if_L1_pos / natural      Ratio  0.93643936
```

**gfoRmula** — user-supplied function passed as an argument:

``` r

treat_if_L1_pos <- function(newdf, pool, intvar, intvals, time_name, t) {
  newdf[, (intvar) := as.integer(newdf[[intvar]] > 0)]
}

gfoRmula::gformula(
  ...,
  intervention1.A = list(treat_if_L1_pos),
  int_descript    = "Treat if L1 > 0"
)
```

Both approaches evaluate the rule within the simulated dataset at each
time step, so the exposure reflects its natural-course draw before the
threshold is applied. The current time step is available inside the
[`dyn_int()`](https://adayim.github.io/causalMed/reference/dyn_int.md)
expression as the tracked `time_var` column
(e.g. `dyn_int(as.numeric(L1 > 0 & time >= 2))`), and history-dependent
rules can be built by maintaining lagged columns through `in_recode` /
`out_recode`. The difference from `gfoRmula` is interface rather than
capability: `gfoRmula` hands a user function the full pooled data
object, whereas `causalMed` expresses the same logic through its
recode-plus-expression pattern.

------------------------------------------------------------------------

## High-Precision Numerical Cross-Validation

With 50,000 Monte Carlo replicates and the same random seed, the two
packages are expected to agree to within 0.0001. The differences shown
below reflect Monte Carlo sampling variability only, not algorithmic
divergence.

``` r

fit_cv_cm <- causalMed::gformula(
  data         = dat, id_var = "id", time_var = "time", base_vars = "V",
  exposure     = "A", models = models_cm,
  intervention = list(natural = NULL, always = 1),
  ref_int      = "natural",
  init_recode  = init_rc, in_recode = in_rc,
  mc_sample    = 50000, R = 1, quiet = TRUE, seed = 20250915
)

fit_cv_gf <- gfoRmula::gformula(
  obs_data = copy(dat), id = "id", time_name = "time",
  time_points = length(unique(dat$time)),
  covnames = c("L1","L2","A"), covtypes = c("normal","binary","binary"),
  basecovs = "V",
  covparams = list(covmodels = c(
    L1 ~ lag1_A + lag1_L1 + V + time,
    L2 ~ lag1_A + lag1_L2 + V + time,
    A  ~ lag1_A + L1 + L2 + V + time
  )),
  ymodel = Y_bin ~ A + L1 + L2, outcome_name = "Y_bin",
  outcome_type = "binary_eof",
  histories = c(gfoRmula::lagged), histvars = list(c("A","L1","L2")),
  intervention1.A = list(gfoRmula::static, rep(1, length(unique(dat$time)))),
  int_descript = "Always treat",
  nsimul = 50000, seed = 20250915
)

cm_vals <- round(fit_cv_cm$effect_size$Est, 5)
gf_vals <- round(fit_cv_gf$result[k == max(fit_cv_gf$result$k)][["g-form mean"]], 5)

data.frame(
  Intervention = c("natural", "always"),
  causalMed    = cm_vals,
  gfoRmula     = gf_vals,
  Difference   = round(abs(cm_vals - gf_vals), 5)
)
#>   Intervention causalMed gfoRmula Difference
#> 1      natural   0.14120  0.14120      0e+00
#> 2       always   0.15073  0.15072      1e-05
```

The differences are well within 0.0001, confirming that `causalMed`’s
g-formula engine is numerically consistent with `gfoRmula` as the
reference implementation.

------------------------------------------------------------------------

## References

1.  Westreich, D., Cole, S. R., Young, J. G., et al. (2012). The
    parametric g-formula to estimate the effect of HAART on incident
    AIDS or death. *Statistics in Medicine*, 31, 2000–2009.

2.  McGrath, S., Lin, V., Zhang, Z., et al. (2020). gfoRmula: An R
    Package for Estimating the Effects of Sustained Treatment Strategies
    via the Parametric g-Formula. *Patterns*, 1, 100008.

3.  Lin, S. H., Young, J. G., Logan, R., & VanderWeele, T. J. (2017).
    Mediation analysis for a survival outcome with time-varying
    exposures, mediators, and confounders. *Statistics in Medicine*,
    36(26), 4153–4166.

4.  Zheng, W., & van der Laan, M. (2017). Longitudinal mediation
    analysis with time-varying mediators and exposures, with application
    to survival outcomes. *Journal of Causal Inference*, 5(2).

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] data.table_1.18.4 causalMed_0.0.4  
#> 
#> loaded via a namespace (and not attached):
#>  [1] future.apply_1.20.2 gtable_0.3.6        jsonlite_2.0.0     
#>  [4] dplyr_1.2.1         compiler_4.6.1      tidyselect_1.2.1   
#>  [7] gfoRmula_1.1.1      stringr_1.6.0       parallel_4.6.1     
#> [10] jquerylib_0.1.4     globals_0.19.1      systemfonts_1.3.2  
#> [13] scales_1.4.0        textshaping_1.0.5   yaml_2.3.12        
#> [16] fastmap_1.2.0       ggplot2_4.0.3       R6_2.6.1           
#> [19] generics_0.1.4      knitr_1.51          htmlwidgets_1.6.4  
#> [22] future_1.70.0       tibble_3.3.1        desc_1.4.3         
#> [25] nnet_7.3-20         bslib_0.11.0        pillar_1.11.1      
#> [28] RColorBrewer_1.1-3  rlang_1.3.0         stringi_1.8.7      
#> [31] cachem_1.1.0        xfun_0.59           fs_2.1.0           
#> [34] sass_0.4.10         S7_0.2.2            otel_0.2.0         
#> [37] cli_3.6.6           progressr_1.0.0     pkgdown_2.2.0      
#> [40] magrittr_2.0.5      digest_0.6.39       grid_4.6.1         
#> [43] lifecycle_1.0.5     vctrs_0.7.3         evaluate_1.0.5     
#> [46] glue_1.8.1          listenv_1.0.0       farver_2.1.2       
#> [49] codetools_0.2-20    ragg_1.5.2          parallelly_1.48.0  
#> [52] rmarkdown_2.31      tools_4.6.1         pkgconfig_2.0.3    
#> [55] htmltools_0.5.9
```
