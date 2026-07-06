
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
  - `mediation_type = "I"`: interventional IDE/IIE (VanderWeele &
    Tchetgen Tchetgen 2017; Lin et al. 2017; Yamamuro et al. 2021) —
    cross-world mediator drawn as a **joint M(1:T) trajectory** by
    row-permuting the reference (a\*) cohort, matching the SAS
    `mGFORMULA` macro and the algorithm of Yamamuro et al. 2021. Does
    not require cross-world independence. **Multiple mediator models**
    are supported under this estimand (Yamamuro 2021 §4); per-mediator
    IIE(M_k) is reported alongside the overall IDE/TE.
  - `mediation_type = "N"`: natural NDE/NIE (Zheng & van der Laan 2017)
    — conditional mediator distribution with exposure swapping; requires
    stronger assumptions (single mediator only)
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
#>   Outcome      : Y_bin  [mean outcome at t = 4, end of follow-up]
#>   Time variable: time  (5 time points: 0 ... 4)
#>   ID variable  : id
#>   Baseline vars: V
#>   Data         : 3,000 individuals, 15,000 observations
#>   MC sample    : 10000
#>   Bootstrap R  : 100
#>   Seed         : 2025
#>   Reference    : natural
#> 
#> --- Mean outcome by intervention --- 
#>    Intervention    Est     Sd 2.5%(pct) 97.5%(pct) 2.5%(norm) 97.5%(norm)
#>          <fctr>  <num>  <num>     <num>      <num>      <num>       <num>
#> 1:      natural 0.1405 0.0067    0.1296     0.1547     0.1275      0.1536
#> 2: always_treat 0.1508 0.0074    0.1380     0.1652     0.1362      0.1654
#> 3:  never_treat 0.0867 0.0139    0.0617     0.1137     0.0596      0.1139
#>   Observed (nonparametric) mean of Y_bin at t = 4 (end of follow-up): 0.1407
#>   (informal model check: compare with the natural-course intervention)
#> 
#> --- Contrasts vs. reference intervention --- 
#>              Intervention  Risk_type Estimate     Sd 2.5%(pct) 97.5%(pct)
#>                    <char>     <char>    <num>  <num>     <num>      <num>
#> 1: always_treat - natural Difference   0.0103 0.0024    0.0054     0.0146
#> 2: always_treat / natural      Ratio   1.0731 0.0170    1.0395     1.1024
#> 3:  never_treat - natural Difference  -0.0538 0.0132   -0.0780    -0.0258
#> 4:  never_treat / natural      Ratio   0.6171 0.0929    0.4514     0.8109
#>    2.5%(norm) 97.5%(norm)
#>         <num>       <num>
#> 1:     0.0055      0.0150
#> 2:     1.0398      1.1065
#> 3:    -0.0797     -0.0279
#> 4:     0.4350      0.7992
#> 
#>   95% CIs: percentile (pct) and normal approximation (norm) from 100 bootstrap replicates.
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
#>   Mediator(s)  : M
#>   Outcome      : Y_bin  [mean outcome at t = 4, end of follow-up]
#>   Time variable: time  (5 time points: 0 ... 4)
#>   ID variable  : id
#>   Baseline vars: V
#>   Data         : 3,000 individuals, 15,000 observations
#>   MC sample    : 10000
#>   Bootstrap R  : 100
#>   n_vw         : 2  (permutation draws averaged per cross-world intervention)
#>   Seed         : 2025
#>   Mediation    : Interventional effects (IDE/IIE) -- Lin et al. (2017)
#> 
#> --- Marginal mean outcome per intervention --- 
#>   Under interventional effects, each intervention draws its mediators from independently-permuted pools (G):
#>   Phi11 = E[Y(a=1, G1)]:  exposure=1, mediators ~ a=1 pool  [reference]
#>   Phi10 = E[Y(a=1, G0)]:  exposure=1, mediators ~ a=0 pool  [cross-world]
#>   Phi00 = E[Y(a=0, G0)]:  exposure=0, mediators ~ a=0 pool  [reference]
#>   nat1/nat0 = E[Y(a=1)]/E[Y(a=0)]:  natural course (used for the total effect)
#>    Intervention    Est     Sd 2.5%(pct) 97.5%(pct) 2.5%(norm) 97.5%(norm)
#>          <char>  <num>  <num>     <num>      <num>      <num>       <num>
#> 1:         nat0 0.0847 0.0135    0.0607     0.1128     0.0582      0.1112
#> 2:         nat1 0.1519 0.0075    0.1383     0.1661     0.1371      0.1667
#> 3:        Phi00 0.0848 0.0135    0.0608     0.1126     0.0584      0.1112
#> 4:        Phi10 0.1456 0.0074    0.1329     0.1586     0.1312      0.1600
#> 5:        Phi11 0.1517 0.0075    0.1389     0.1662     0.1370      0.1664
#>   Observed (nonparametric) mean of Y_bin at t = 4 (end of follow-up): 0.1407
#>   (informal benchmark; interventions fix the exposure, so exact agreement is not expected)
#> 
#> --- Effect decomposition --- 
#>   Direct effect (IDE)   = Phi10 - Phi00
#>   Indirect effect (IIE) = Phi11 - Phi10   (sequential per mediator when N>=2)
#>   IDE + IIE             = Phi11 - Phi00    (interventional overall effect)
#>   Total effect (TE)     = nat1 - nat0      (natural plug-in g-formula)
#>   TE - (Direct+Indirect)= mediated-interaction residual (TE - overall)
#>   Mediation Prop.       = (Total - Direct) / Total  (percentage; RR not applicable)
#>   RD = risk difference;  RR = risk ratio
#>                                   Effect     RD     RR Sd(RD) RD 2.5%(pct)
#>                                   <char>  <num>  <num>  <num>        <num>
#> 1:                       Indirect effect 0.0061 1.0421 0.0021       0.0025
#> 2:                         Direct effect 0.0608 1.7168 0.0150       0.0287
#> 3:                          Total effect 0.0672 1.7931 0.0153       0.0323
#> 4:              TE - (Direct + Indirect) 0.0003     NA 0.0004      -0.0009
#> 5:                  Mediation Proportion 9.5299     NA 3.7089       4.0617
#> 6: Mediation Proportion (multiplicative) 9.1585     NA 3.7171       4.0309
#>    RD 97.5%(pct) Sd(RR) RR 2.5%(pct) RR 97.5%(pct) RD 2.5%(norm) RD 97.5%(norm)
#>            <num>  <num>        <num>         <num>         <num>          <num>
#> 1:        0.0108 0.0149       1.0178        1.0732        0.0020         0.0103
#> 2:        0.0856 0.2924       1.2632        2.3606        0.0313         0.0903
#> 3:        0.0916 0.3076       1.2961        2.4283        0.0371         0.0973
#> 4:        0.0007     NA           NA            NA       -0.0005         0.0011
#> 5:       17.6448     NA           NA            NA        2.2606        16.7991
#> 6:       17.7527     NA           NA            NA        1.8732        16.4438
#>    RR 2.5%(norm) RR 97.5%(norm)
#>            <num>          <num>
#> 1:        1.0129         1.0713
#> 2:        1.1437         2.2899
#> 3:        1.1903         2.3959
#> 4:            NA             NA
#> 5:            NA             NA
#> 6:            NA             NA
#> 
#>   95% CIs: percentile (pct) and normal approximation (norm) from 100 bootstrap replicates.
```

The `estimate` component of `fit_med` contains:

| Effect | Definition |
|----|----|
| Indirect effect | Q(1,1) − Q(1,0): effect through the mediator pathway |
| Direct effect | Q(1,0) − Q(0,0): effect not through the mediator |
| Total effect | natural plug-in g-formula contrast E\[Y₁\] − E\[Y₀\] |
| TE − (Direct + Indirect) | mediated-interaction residual (interventional only; exactly 0, and omitted, for natural effects) |
| Mediation Proportion | (Total − Direct) / Total × 100% (additive) |
| Mediation Proportion (multiplicative) | RR-scale proportion mediated (Lin et al. 2017, Table 2) |

With `R > 1`, `estimate` additionally carries bootstrap SEs and
percentile/normal CIs on the RD and RR scales, and the returned object
gains `boot_estimates` (the per-replicate draws), `data_summary`, and an
`observed` nonparametric benchmark that `print()` displays next to the
simulated means.

For `mediation_type = "I"`, Q(a₁, a₂) is the mean outcome when exposure
is set to a₁ and each mediator is a **stochastic draw** from its
marginal distribution under exposure a₂. The direct and indirect effects
then sum to the *interventional overall effect* Q(1,1) − Q(0,0), which
generally differs from the total effect; their gap is the residual row.
For `mediation_type = "N"` the references use each subject’s own
mediator and the decomposition sums exactly to the total effect (no
residual row).

Under `mediation_type = "I"`, the mediators in **every** intervention —
the references Q(0,0) and Q(1,1) as well as the cross-world Q(1,0) — are
drawn as **joint trajectories**: a natural-course cohort is simulated
under each a\*, and the intervention assigns each subject the full
M(1:T) of a randomly permuted reference-cohort individual (each mediator
permuted independently). This matches the SAS `mGFORMULA` macro (Lin et
al. 2017, eAppendix) and Yamamuro et al. 2021 (Figure 3, step 3).
Mediator values are not survival-weighted; the full reference cohort is
used. See `?mediation` for details.

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
Variables used in models as lags (e.g., `lag1_A`) must be pre-computed
and present in the data, or created via `init_recode` / `in_recode`
hooks. The `nonsurvivaldata` bundled with the package illustrates the
required structure.

``` r
data("nonsurvivaldata", package = "causalMed")
head(nonsurvivaldata)
#>       id  time         V         L1    L2     A          M    Y_cont Y_bin
#>    <num> <num>     <num>      <num> <num> <num>      <num>     <num> <num>
#> 1:     1     0 0.5218287  0.5527842     0     1 -0.2196759        NA    NA
#> 2:     1     1 0.5218287 -0.1758057     1     1  0.4719272        NA    NA
#> 3:     1     2 0.5218287  0.5903320     0     1 -0.2237722        NA    NA
#> 4:     1     3 0.5218287  1.3174098     0     1  0.3878269        NA    NA
#> 5:     1     4 0.5218287  0.7974775     0     1 -0.4163290 0.1348551     0
#> 6:     2     0 0.4066671  0.5492001     1     1  1.0580998        NA    NA
#>    lag1_A    lag1_L1 lag1_L2     lag1_M
#>     <num>      <num>   <num>      <num>
#> 1:     NA         NA      NA         NA
#> 2:      1  0.5527842       0 -0.2196759
#> 3:      1 -0.1758057       1  0.4719272
#> 4:      1  0.5903320       0 -0.2237722
#> 5:      1  1.3174098       0  0.3878269
#> 6:     NA         NA      NA         NA
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

5.  VanderWeele, T. J., & Tchetgen Tchetgen, E. J. (2017). Mediation
    analysis with time varying exposures and mediators. *Journal of the
    Royal Statistical Society: Series B*, 79(3), 917–938.
    [doi:10.1111/rssb.12194](https://doi.org/10.1111/rssb.12194)

6.  Yamamuro, S., Shinozaki, T., Iimuro, S., & Matsuyama, Y. (2021).
    Mediational g-formula for time-varying treatment and
    repeated-measured multiple mediators: Application to atorvastatin’s
    effect on cardiovascular disease via cholesterol lowering and
    anti-inflammatory actions in elderly type 2 diabetics. *Statistical
    Methods in Medical Research*, 30(8), 1782–1799.
    [doi:10.1177/09622802211025988](https://doi.org/10.1177/09622802211025988)
