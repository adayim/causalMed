
<!-- README.md is generated from README.Rmd. Please edit that file -->

# causalMed

**Causal Mediation Analysis with Time-Varying Exposures, Mediators, and
Confounders**

`causalMed` implements the **parametric g-formula** for total effect
estimation and extends it with the **survival mediational g-formula** to
decompose causal effects into direct and indirect components. The
package supports:

- **Interventional direct and indirect effects** (IDE/IIE) — Lin et
  al. (2017)
- **Natural direct and indirect effects** (NDE/NIE) — Zheng & van der
  Laan (2017)

Both approaches handle time-varying exposures, mediators, and
confounders in longitudinal data, including survival outcomes. The
standard g-formula component (total effect estimation) is
cross-validated against the `gfoRmula` CRAN package.

## Installation

``` r
# Install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("adayim/causalMed")
```

## Key Features

- **G-formula for total effects**: Monte Carlo simulation under
  user-defined static or dynamic (threshold) interventions for binary,
  continuous, and survival outcomes
- **Mediational g-formula**: decomposes total effects into direct and
  indirect pathways with time-varying mediators and exposure-induced
  mediator–outcome confounders
- **Two mediation estimands**:
  - `mediation_type = "I"`: interventional IDE/IIE (Lin et al. 2017) —
    marginal mediator distribution via random permutation; does not
    require cross-world independence
  - `mediation_type = "N"`: natural NDE/NIE (Zheng & van der Laan 2017)
    — conditional mediator distribution with exposure swapping; requires
    stronger assumptions
- **Flexible model specification**: logistic regression (binary), linear
  regression (normal), multinomial logistic (categorical), and custom
  simulation functions
- **Recode hooks** (`init_recode`, `in_recode`, `out_recode`) for lag
  creation, cumulative variables, and other within-loop transformations
- **Bootstrap confidence intervals**: non-parametric bootstrap
  resampling individuals (preserving longitudinal structure), with
  percentile and normal-approximation CIs
- **Parallel bootstrap**: plug in `future::plan(multisession)` on
  Windows or `future::plan(multicore)` on Unix to parallelise across
  bootstrap replicates

## Usage

### Total effect: standard parametric g-formula

``` r
library(causalMed)

data("nonsurvivaldata", package = "causalMed")

# Manage lagged variables
init_rc <- recodes(lag1_A  = 0,   # At t=0, all lags initialised to 0
                   lag1_L1 = 0,
                   lag1_L2 = 0)

in_rc   <- recodes(lag1_A  = A,   # At each subsequent step, copy current values
                   lag1_L1 = L1,
                   lag1_L2 = L2)

# Specify models in temporal order: L1 → L2 → A → Y
m_L1 <- spec_model(L1    ~ lag1_A + lag1_L1 + V + time,
                   var_type = "normal",  mod_type = "covariate")
m_L2 <- spec_model(L2    ~ lag1_A + lag1_L2 + V + time,
                   var_type = "binary",  mod_type = "covariate")
m_A  <- spec_model(A     ~ lag1_A + L1 + L2 + V + time,
                   var_type = "binary",  mod_type = "exposure")
m_Y  <- spec_model(Y_bin ~ A + L1 + L2,
                   var_type = "binary",  mod_type = "outcome")

models_bin <- list(m_L1, m_L2, m_A, m_Y)

# Define intervention strategies
# NULL  = natural course (draw exposure from its fitted model)
# 1 / 0 = always treat / never treat
ints <- list(natural = NULL, always_treat = 1, never_treat = 0)

# ── 3. Run g-formula ────────────────────────────────────────────────────────
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
  R           = 100,  # set R > 1 for bootstrap CIs; kept low here for speed
  quiet       = TRUE,
  seed        = 2025
)
#> Warning: package 'future' was built under R version 4.5.3

print(fit_bin)
#> Call:
#> gformula(data = nonsurvivaldata, id_var = "id", base_vars = "V", 
#>     exposure = "A", time_var = "time", models = models_bin, intervention = ints, 
#>     ref_int = "natural", init_recode = init_rc, in_recode = in_rc, 
#>     mc_sample = 10000, R = 100, quiet = TRUE, seed = 2025)
#> 
#> --- Analysis setup ---
#>   Exposure     : A
#>   Time variable: time
#>   ID variable  : id
#>   Baseline vars: V
#>   MC sample    : 10000
#>   Bootstrap R  : 100
#>   Seed         : 2025
#>   Reference    : natural
#> 
#> --- Mean outcome by intervention --- 
#>    Intervention    Est     Sd RD 2.5%(pct) RD 97.5%(pct) RD 2.5%(norm)
#>          <fctr>  <num>  <num>        <num>         <num>         <num>
#> 1:      natural 0.1405 0.0067       0.1296        0.1547        0.1275
#> 2: always_treat 0.1508 0.0074       0.1380        0.1652        0.1362
#> 3:  never_treat 0.0867 0.0139       0.0617        0.1137        0.0596
#>    RD 97.5%(norm)
#>             <num>
#> 1:         0.1536
#> 2:         0.1654
#> 3:         0.1139
#> 
#> --- Contrasts vs. reference intervention --- 
#>              Intervention  Risk_type Estimate     Sd RD 2.5%(pct) RD 97.5%(pct)
#>                    <char>     <char>    <num>  <num>        <num>         <num>
#> 1: always_treat - natural Difference   0.0103 0.0024       0.0054        0.0146
#> 2: always_treat / natural      Ratio   1.0731 0.0170       1.0395        1.1024
#> 3:  never_treat - natural Difference  -0.0538 0.0132      -0.0780       -0.0258
#> 4:  never_treat / natural      Ratio   0.6171 0.0929       0.4514        0.8109
#>    RD 2.5%(norm) RD 97.5%(norm)
#>            <num>          <num>
#> 1:        0.0055         0.0150
#> 2:        1.0398         1.1065
#> 3:       -0.0797        -0.0279
#> 4:        0.4350         0.7992
```

### Mediation analysis: interventional IDE/IIE (Lin et al. 2017)

``` r
library(causalMed)

data("nonsurvivaldata", package = "causalMed")

# Model list must include a mediator model (mod_type = "mediator")
# List order must reflect the temporal DAG: A -> L -> M -> Y
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

fit_med <- mediation(
  data           = nonsurvivaldata,
  id_var         = "id",
  time_var       = "time",
  base_vars      = "V",
  exposure       = "A",
  outcome        = "Y_bin",
  models         = models_med,
  init_recode    = init_med,
  in_recode      = in_med,
  mediation_type = "I",     # Interventional IDE/IIE
  mc_sample      = 10000,
  R              = 100,
  quiet          = TRUE,
  seed           = 2025
)

print(fit_med)
#> Call:
#> mediation(data = nonsurvivaldata, id_var = "id", base_vars = "V", 
#>     exposure = "A", outcome = "Y_bin", time_var = "time", models = models_med, 
#>     init_recode = init_med, in_recode = in_med, mc_sample = 10000, 
#>     mediation_type = "I", R = 100, quiet = TRUE, seed = 2025)
#> 
#> --- Analysis setup ---
#>   Exposure     : A
#>   Time variable: time
#>   ID variable  : id
#>   Baseline vars: V
#>   MC sample    : 10000
#>   Bootstrap R  : 100
#>   Seed         : 2025
#>   Mediation    : Interventional effects (IDE/IIE) -- Lin et al. (2017)
#> 
#> --- Marginal mean outcome per arm (Q-functionals) --- 
#>   Phi11 = E[Y(a=1, M(1))]:  exposure=1, mediator under a=1
#>   Phi10 = E[Y(a=1, M(0))]:  exposure=1, mediator under a=0  [cross-world]
#>   Phi00 = E[Y(a=0, M(0))]:  exposure=0, mediator under a=0
#>    Intervention    Est     Sd RD 2.5%(pct) RD 97.5%(pct) RD 2.5%(norm)
#>          <fctr>  <num>  <num>        <num>         <num>         <num>
#> 1:        Phi11 0.1519 0.0075       0.1383        0.1661        0.1371
#> 2:        Phi00 0.0847 0.0135       0.0607        0.1128        0.0582
#> 3:        Phi10 0.1460 0.0073       0.1329        0.1589        0.1317
#>    RD 97.5%(norm)
#>             <num>
#> 1:         0.1667
#> 2:         0.1112
#> 3:         0.1604
#> 
#> --- Effect decomposition --- 
#>   Total effect    = Phi11 - Phi00 =  E[Y(1,M(1))] - E[Y(0,M(0))]
#>   Direct effect   = Phi10 - Phi00 =  E[Y(1,M(0))] - E[Y(0,M(0))]
#>   Indirect effect = Phi11 - Phi10 =  E[Y(1,M(1))] - E[Y(1,M(0))]
#>   Mediation Prop. = Indirect / Total  (as a percentage; RR not applicable)
#>   RD = risk difference;  RR = risk ratio
#>                  Effect     RD     RR     Sd RD 2.5%(pct) RD 97.5%(pct)  Sd_RR
#>                  <char>  <num>  <num>  <num>        <num>         <num>  <num>
#> 1:      Indirect effect 0.0059 1.0401 0.0021       0.0024        0.0105 0.0148
#> 2:        Direct effect 0.0613 1.7239 0.0151       0.0287        0.0854 0.2938
#> 3:         Total effect 0.0672 1.7931 0.0153       0.0323        0.0916 0.3076
#> 4: Mediation Proportion 8.7231     NA 3.6941       4.3050       17.5163     NA
#>    RR 2.5%(pct) RR 97.5%(pct) RD 2.5%(norm) RD 97.5%(norm) RR 2.5%(norm)
#>           <num>         <num>         <num>          <num>         <num>
#> 1:       1.0164        1.0729        0.0018         0.0099        1.0112
#> 2:       1.2612        2.3621        0.0318         0.0909        1.1480
#> 3:       1.2961        2.4283        0.0371         0.0973        1.1903
#> 4:           NA            NA        1.4829        15.9634            NA
#>    RR 97.5%(norm)
#>             <num>
#> 1:         1.0690
#> 2:         2.2997
#> 3:         2.3959
#> 4:             NA
```

The `estimate` component of `fit_med` contains:

| Effect               | Definition                                           |
|----------------------|------------------------------------------------------|
| Indirect effect      | Q(1,1) − Q(1,0): effect through the mediator pathway |
| Direct effect        | Q(1,0) − Q(0,0): effect not through the mediator     |
| Total effect         | Q(1,1) − Q(0,0): sum of direct and indirect          |
| Mediation proportion | Indirect / Total × 100%                              |

where Q(a₁, a₂) is the mean outcome when exposure is set to a₁ and the
mediator distribution is drawn from the distribution it would have had
under exposure a₂.

### Natural NDE/NIE (Zheng & van der Laan 2017)

Change one argument:

``` r
fit_natural <- mediation(
  ...,                     # same arguments as above
  mediation_type = "N"    # natural NDE/NIE (Zheng & van der Laan 2017)
)
```

Natural effects condition the mediator model on the individual’s own
covariate history but evaluate it at the alternative exposure level
(exposure swapping). They require stronger sequential
no-unmeasured-confounding assumptions than interventional effects.

### Enabling parallel bootstrap

``` r
library(future)
plan(multisession)   # Windows; use plan(multicore) on Unix/macOS

fit_parallel <- mediation(..., R = 500)

plan(sequential)     # reset after use
```

## Data format

Input data must be in **long format** with one row per subject per time
point. The time variable should be an ordered integer starting at 0.
Variables used in models as lags (e.g., `A_lag1`) must be pre-computed
and present in the data, or created via `init_recode` / `in_recode`
hooks. The `nonsurvivaldata` bundled with the package illustrates the
required structure.

``` r
data("nonsurvivaldata", package = "causalMed")
head(nonsurvivaldata)
#>   id time V          L A         M Y M_lag1    L_lag1 A_lag1
#> 1  1    0 1 -0.4557851 0 0.1847927 0      0  0.000000      0
#> 2  1    1 1  0.4651968 0 0.9559358 0      0 -0.455785      0
#> 3  1    2 1  0.2498987 1 1.3503895 0      1  0.465197      0
```

## References

1.  Westreich, D., Cole, S. R., Young, J. G., et al. (2012). The
    parametric g-formula to estimate the effect of HAART on incident
    AIDS or death. *Statistics in Medicine*, 31, 2000–2009.
    [doi:10.1002/sim.5316](https://doi.org/10.1002/sim.5316)

2.  McGrath, S., Lin, V., Zhang, Z., et al. (2020). gfoRmula: An R
    Package for Estimating the Effects of Sustained Treatment Strategies
    via the Parametric g-Formula. *Patterns*, 1, 100008.
    [doi:10.1016/j.patter.2020.100008](https://doi.org/10.1016/j.patter.2020.100008)

3.  Lin, S. H., Young, J. G., Logan, R., & VanderWeele, T. J. (2017).
    Mediation analysis for a survival outcome with time-varying
    exposures, mediators, and confounders. *Statistics in Medicine*,
    36(26), 4153–4166.
    [doi:10.1002/sim.7426](https://doi.org/10.1002/sim.7426)

4.  Zheng, W., & van der Laan, M. (2017). Longitudinal mediation
    analysis with time-varying mediators and exposures, with application
    to survival outcomes. *Journal of Causal Inference*, 5(2).
    [doi:10.1515/jci-2016-0006](https://doi.org/10.1515/jci-2016-0006)
