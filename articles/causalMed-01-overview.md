# Getting Started with causalMed

## Introduction

**causalMed** applies the parametric g-formula (Robins 1986) to
longitudinal data with time-varying exposures, mediators and
confounders, including confounders affected by earlier exposure. It has
two main functions:

- [`gformula()`](https://adayim.github.io/causalMed/reference/gformula.md)
  estimates **total effects**: the mean outcome or risk had everyone
  followed a given exposure strategy (Westreich et al. 2012; McGrath et
  al. 2020).
- [`mediation()`](https://adayim.github.io/causalMed/reference/mediation.md)
  splits the effect of an exposure into **direct and indirect effects**
  through one or more mediators, either interventional (Lin et al.
  2017. or natural (Zheng & van der Laan 2017).

This vignette covers what every analysis needs: the data format, model
specification, lagged variables, and a quick start for each function.
The other vignettes go further:

| Vignette | Contents |
|----|----|
| [`vignette("causalMed-02-mediation")`](https://adayim.github.io/causalMed/articles/causalMed-02-mediation.md) | The two mediation estimands, survival outcomes, multiple mediators, censoring, natural effects. |
| [`vignette("causalMed-03-gformula")`](https://adayim.github.io/causalMed/articles/causalMed-03-gformula.md) | Total effects, dynamic interventions, custom distributions, bootstrap intervals, a published replication. |
| [`vignette("causalMed-04-vs-gfoRmula")`](https://adayim.github.io/causalMed/articles/causalMed-04-vs-gfoRmula.md) | Comparison with the CRAN package `gfoRmula`. |

``` r

library(causalMed)
library(data.table)
```

## Data

Data must be in **long format**: one row per subject per time point,
with a numeric time variable.

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

`nonsurvivaldata` follows 3,000 subjects at times 0 to 4:

| Variable   | Role                                   |
|------------|----------------------------------------|
| `id`       | Subject identifier                     |
| `time`     | Time index (0 to 4)                    |
| `V`        | Baseline covariate                     |
| `A`        | Binary exposure                        |
| `L1`, `L2` | Continuous and binary confounders      |
| `M`        | Continuous mediator                    |
| `Y_bin`    | Binary outcome at the end of follow-up |

Within each time point the exposure comes first, the confounders respond
to it, then the mediator, then the outcome (see
[`?nonsurvivaldata`](https://adayim.github.io/causalMed/reference/nonsurvivaldata.md)).

## Specifying models

Each time-varying variable needs a model, created with
[`spec_model()`](https://adayim.github.io/causalMed/reference/spec_model.md).
The models are collected in a list **in the order the variables are
generated** within a time point.

``` r

spec_model(
  formula,           # response ~ predictors
  var_type,          # how values are simulated
  mod_type,          # role in the causal structure
  subset = NULL,     # optional: restrict the model to some rows
  recode = NULL      # optional: recodes() applied before this model
)
```

| `var_type` | Model and simulated values |
|----|----|
| `"normal"` | Linear regression; Gaussian draws, clipped to the observed range by default (`truncate = TRUE`) |
| `"binary"` | Logistic regression; Bernoulli draws |
| `"categorical"` | Multinomial logistic regression ([`nnet::multinom`](https://rdrr.io/pkg/nnet/man/multinom.html)) |
| `"custom"` | Your own fitting and simulation functions (see [`vignette("causalMed-03-gformula")`](https://adayim.github.io/causalMed/articles/causalMed-03-gformula.md)) |

| `mod_type` | Role |
|----|----|
| `"exposure"` | The variable intervened on |
| `"covariate"` | Time-varying confounder |
| `"mediator"` | Mediator (needed by [`mediation()`](https://adayim.github.io/causalMed/reference/mediation.md)) |
| `"outcome"` | Outcome at the end of follow-up |
| `"survival"` | Discrete-time event indicator |
| `"censor"` | Loss-to-follow-up indicator |

A confounder affected by the exposure at the same time point must come
after the exposure model and before the outcome model.

## Lagged and derived variables

The simulation starts from the subject identifier and the baseline
covariates only. Every lag, counter or other derived variable a model
uses must therefore be created with
[`recodes()`](https://adayim.github.io/causalMed/reference/recodes.md),
even if it already exists in the data; otherwise the run stops with an
error such as `object 'lag1_A' not found`.

| Hook | When it runs | Typical use |
|----|----|----|
| `init_recode` | First time step, before the models | Initial values of lags |
| `in_recode` | Start of each later time step, before the models | Update lags |
| `out_recode` | End of each time step, after the models | Cumulative counts, absorbing states |

``` r

init_rc <- recodes(lag1_A = 0, lag1_L1 = 0, lag1_L2 = 0)   # first time step
in_rc   <- recodes(lag1_A = A, lag1_L1 = L1, lag1_L2 = L2) # later steps
```

Lag columns do not belong in `base_vars`: they are not time-fixed, and a
warning is given if they are listed there.

## Quick start: total effects with `gformula()`

Models in the order **A → L1 → L2 → Y**, and three strategies: the
natural course, always exposed, and never exposed.

``` r

models_te <- list(
  spec_model(A     ~ V + lag1_A + lag1_L1 + lag1_L2 + time,
             var_type = "binary", mod_type = "exposure"),
  spec_model(L1    ~ V + A + lag1_L1 + time,
             var_type = "normal", mod_type = "covariate"),
  spec_model(L2    ~ V + A + lag1_L2 + time,
             var_type = "binary", mod_type = "covariate"),
  spec_model(Y_bin ~ V + A + L1 + L2,
             var_type = "binary", mod_type = "outcome")
)

fit_te <- gformula(
  data         = nonsurvivaldata,
  id_var       = "id",
  time_var     = "time",
  base_vars    = "V",
  exposure     = "A",
  models       = models_te,
  intervention = list(natural = NULL, always_treat = 1, never_treat = 0),
  ref_int      = "natural",
  init_recode  = init_rc,
  in_recode    = in_rc,
  mc_sample    = 10000,
  R            = 1,        # set R > 1 for bootstrap confidence intervals
  quiet        = TRUE,
  seed         = 2025
)

fit_te$effect_size   # mean outcome under each strategy
#>    Intervention       Est
#>          <fctr>     <num>
#> 1:      natural 0.2349247
#> 2: always_treat 0.2522388
#> 3:  never_treat 0.1066566
fit_te$estimate      # risk difference and ratio against the natural course
#>              Intervention  Risk_type    Estimate
#>                    <char>     <char>       <num>
#> 1: always_treat - natural Difference  0.01731407
#> 2: always_treat / natural      Ratio  1.07370051
#> 3:  never_treat - natural Difference -0.12826810
#> 4:  never_treat / natural      Ratio  0.45400334
```

## Quick start: mediation with `mediation()`

Add a mediator model and name the outcome. The default estimand is the
interventional direct and indirect effect (Lin et al. 2017):

``` r

init_med <- recodes(lag1_A = 0, lag1_L1 = 0, lag1_L2 = 0, lag1_M = 0)
in_med   <- recodes(lag1_A = A, lag1_L1 = L1, lag1_L2 = L2, lag1_M = M)

models_med <- list(
  spec_model(A     ~ V + lag1_A + lag1_L1 + lag1_L2 + time,
             var_type = "binary", mod_type = "exposure"),
  spec_model(L1    ~ V + A + lag1_L1 + time,
             var_type = "normal", mod_type = "covariate"),
  spec_model(L2    ~ V + A + lag1_L2 + time,
             var_type = "binary", mod_type = "covariate"),
  spec_model(M     ~ V + A + L1 + L2 + lag1_M + time,
             var_type = "normal", mod_type = "mediator"),
  spec_model(Y_bin ~ V + A + M + L1 + L2,
             var_type = "binary", mod_type = "outcome")
)

fit_med <- mediation(
  data        = nonsurvivaldata,
  id_var      = "id",
  time_var    = "time",
  base_vars   = "V",
  exposure    = "A",
  outcome     = "Y_bin",
  models      = models_med,
  init_recode = init_med,
  in_recode   = in_med,
  mc_sample   = 10000,
  R           = 1,
  quiet       = TRUE,
  seed        = 2025
)

fit_med$estimate
#>                      Effect           RD       RR
#>                      <char>        <num>    <num>
#> 1:          Indirect effect  0.068268280 1.413236
#> 2:            Direct effect  0.087924817 2.137753
#> 3:             Total effect  0.161254808 2.678705
#> 4: TE - (Direct + Indirect)  0.005061711       NA
#> 5:     Mediation Proportion 43.707616415       NA
```

The table gives the direct effect, the indirect effect through `M`, the
total effect, the decomposition residual (the direct and indirect
effects sum to the interventional overall effect, not to the total
effect) and the proportion mediated.
[`vignette("causalMed-02-mediation")`](https://adayim.github.io/causalMed/articles/causalMed-02-mediation.md)
explains each row, the natural effects of `mediation_type = "N"`, and
how to choose between the two.

## Parallel bootstrap

With `R > 1`, bootstrap replicates run through `future.apply`. Set a
parallel plan before the call:

``` r

future::plan(future::multisession)
fit <- gformula(..., R = 500)
future::plan(future::sequential)
```

## Causal assumptions

The g-formula identifies these effects under assumptions that the
package cannot check (Robins 1986; Westreich et al. 2012; Keil et
al. 2014):

1.  **Consistency**: the observed outcome equals the potential outcome
    under the observed exposure history.
2.  **Positivity**: each exposure level has a non-zero probability for
    every covariate history that occurs under the intervention.
3.  **Sequential exchangeability**: no unmeasured confounding at each
    time point, given the measured past.

The natural effects of `mediation_type = "N"` are identified under
sequential randomization and positivity (Zheng & van der Laan 2017,
Lemma 1). Reading them as **individual-level** natural effects
additionally requires a cross-world assumption that is not expected to
hold when a mediator-outcome confounder is affected by the exposure, as
`L1` and `L2` are here (Avin, Shpitser & Pearl 2005; VanderWeele &
Tchetgen Tchetgen 2017). The interventional effects of
`mediation_type = "I"` do not need it.

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
- McGrath, S., Lin, V., Zhang, Z., et al. (2020). gfoRmula: An R package
  for estimating the effects of sustained treatment strategies via the
  parametric g-formula. *Patterns*, 1, 100008.
- Robins, J. M. (1986). A new approach to causal inference in mortality
  studies with a sustained exposure period: application to control of
  the healthy worker survivor effect. *Mathematical Modelling*, 7(9–12),
  1393–1512.
- VanderWeele, T. J., & Tchetgen Tchetgen, E. J. (2017). Mediation
  analysis with time varying exposures and mediators. *Journal of the
  Royal Statistical Society: Series B*, 79(3), 917–938.
- Westreich, D., Cole, S. R., Young, J. G., et al. (2012). The
  parametric g-formula to estimate the effect of highly active
  antiretroviral therapy on incident AIDS or death. *Statistics in
  Medicine*, 31, 2000–2009.
- Zheng, W., & van der Laan, M. (2017). Longitudinal mediation analysis
  with time-varying mediators and exposures, with application to
  survival outcomes. *Journal of Causal Inference*, 5(2), 20160006.
