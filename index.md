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
    IIE(M_k) is reported alongside the overall IDE/TE. On what a
    non-zero interventional indirect effect does and does not establish,
    see Miles (2023).
  - `mediation_type = "N"`: natural NDE/NIE (Zheng & van der Laan 2017)
    — conditional mediator distribution with exposure swapping; requires
    stronger assumptions (single mediator only). Natural effects are
    **not identifiable** when a confounder of the mediator–outcome
    relationship is itself affected by the exposure (Avin, Shpitser &
    Pearl 2005; VanderWeele & Tchetgen Tchetgen 2017);
    [`mediation()`](https://adayim.github.io/causalMed/reference/mediation.md)
    detects this, warns, and repeats the caveat in
    [`print()`](https://rdrr.io/r/base/print.html). VanderWeele &
    Tchetgen Tchetgen (2017) propose the interventional estimand for
    that setting.
- **Two estimators for natural effects**: the parametric g-formula
  plug-in (`estimator = "gcomp"`, the default) or a **targeted maximum
  likelihood estimator** (`estimator = "tmle"`, Zheng & van der Laan
  2017 §4.3), for which that paper establishes multiple robustness —
  consistency when only certain subsets of the nuisance models are
  correct — with Wald confidence intervals from the efficient influence
  curve, so no bootstrap is needed. Available for `mediation_type = "N"`
  only.
- **Flexible model specification**: logistic regression (binary), linear
  regression (normal), multinomial logistic (categorical), and custom
  simulation functions. Numeric draws are clipped to the observed range
  of the response by default; use `spec_model(truncate = FALSE)` to draw
  from the untruncated fitted distribution.
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

# Specify models in temporal order: A → L1 → L2 → Y (exposure first,
# current A in the confounder models — matching the documented DGP)
m_A  <- spec_model(A     ~ V + lag1_A + lag1_L1 + lag1_L2 + time,
                   var_type = "binary",  mod_type = "exposure")
m_L1 <- spec_model(L1    ~ V + A + lag1_L1 + time,
                   var_type = "normal",  mod_type = "covariate")
m_L2 <- spec_model(L2    ~ V + A + lag1_L2 + time,
                   var_type = "binary",  mod_type = "covariate")
m_Y  <- spec_model(Y_bin ~ V + A + L1 + L2,
                   var_type = "binary",  mod_type = "outcome")

models_bin <- list(m_A, m_L1, m_L2, m_Y)

# Define intervention strategies
# NULL  = natural course (draw exposure from its fitted model)
# 1 / 0 = always treat / never treat
ints <- list(natural = NULL, always_treat = 1, never_treat = 0)

# Run g-formula
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
#> 1:      natural 0.2349 0.0081    0.2207     0.2487     0.2191      0.2508
#> 2: always_treat 0.2522 0.0088    0.2375     0.2681     0.2349      0.2696
#> 3:  never_treat 0.1067 0.0163    0.0790     0.1394     0.0748      0.1385
#>   Observed (nonparametric) mean of Y_bin at t = 4 (end of follow-up): 0.2333
#>   (informal model check: compare with the natural-course intervention)
#> 
#> --- Contrasts vs. reference intervention --- 
#>              Intervention  Risk_type Estimate     Sd 2.5%(pct) 97.5%(pct)
#>                    <char>     <char>    <num>  <num>     <num>      <num>
#> 1: always_treat - natural Difference   0.0173 0.0023    0.0125     0.0212
#> 2: always_treat / natural      Ratio   1.0737 0.0097    1.0559     1.0922
#> 3:  never_treat - natural Difference  -0.1283 0.0164   -0.1549    -0.0908
#> 4:  never_treat / natural      Ratio   0.4540 0.0677    0.3431     0.5975
#>    2.5%(norm) 97.5%(norm)
#>         <num>       <num>
#> 1:     0.0128      0.0219
#> 2:     1.0546      1.0928
#> 3:    -0.1604     -0.0961
#> 4:     0.3213      0.5867
#> 
#>   95% CIs: percentile (pct) and normal approximation (norm) from 100 bootstrap replicates.
```

The list order must match your assumed data-generating process — here
the exposure is decided first within each period and the confounders
respond to it, as documented in
[`?nonsurvivaldata`](https://adayim.github.io/causalMed/reference/nonsurvivaldata.md).

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
#> 1:         nat0 0.0961 0.0151    0.0695     0.1262     0.0664      0.1257
#> 2:         nat1 0.2573 0.0090    0.2400     0.2719     0.2396      0.2750
#> 3:        Phi00 0.0773 0.0132    0.0547     0.1046     0.0514      0.1032
#> 4:        Phi10 0.1652 0.0102    0.1501     0.1883     0.1453      0.1851
#> 5:        Phi11 0.2335 0.0098    0.2157     0.2503     0.2143      0.2527
#>   Observed (nonparametric) mean of Y_bin at t = 4 (end of follow-up): 0.2333
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
#>                                   Effect      RD     RR Sd(RD) RD 2.5%(pct)
#>                                   <char>   <num>  <num>  <num>        <num>
#> 1:                       Indirect effect  0.0683 1.4132 0.0085       0.0521
#> 2:                         Direct effect  0.0879 2.1378 0.0159       0.0571
#> 3:                          Total effect  0.1613 2.6787 0.0179       0.1190
#> 4:              TE - (Direct + Indirect)  0.0051     NA 0.0018      -0.0001
#> 5:                  Mediation Proportion 45.4746     NA 5.6305      34.8497
#> 6: Mediation Proportion (multiplicative) 43.7076     NA 5.9414      34.0057
#>    RD 97.5%(pct) Sd(RR) RR 2.5%(pct) RR 97.5%(pct) RD 2.5%(norm) RD 97.5%(norm)
#>            <num>  <num>        <num>         <num>         <num>          <num>
#> 1:        0.0819 0.0652       1.2857        1.5185        0.0516         0.0849
#> 2:        0.1191 0.3897       1.5371        2.9167        0.0568         0.1190
#> 3:        0.1872 0.4424       1.9739        3.6034        0.1261         0.1964
#> 4:        0.0066     NA           NA            NA        0.0015         0.0086
#> 5:       56.0160     NA           NA            NA       34.4391        56.5101
#> 6:       55.8100     NA           NA            NA       32.0626        55.3526
#>    RR 2.5%(norm) RR 97.5%(norm)
#>            <num>          <num>
#> 1:        1.2854         1.5411
#> 2:        1.3739         2.9016
#> 3:        1.8117         3.5457
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
`observed` nonparametric benchmark that
[`print()`](https://rdrr.io/r/base/print.html) displays next to the
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
used. See
[`?mediation`](https://adayim.github.io/causalMed/reference/mediation.md)
for details.

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
no-unmeasured-confounding assumptions than interventional effects. In
particular, natural direct and indirect effects are **not identifiable**
when a mediator–outcome confounder is itself affected by prior exposure
(Avin, Shpitser & Pearl 2005; VanderWeele & Tchetgen Tchetgen 2017).
[`mediation()`](https://adayim.github.io/causalMed/reference/mediation.md)
scans the model formulas for covariates modelled as exposure-affected,
warns at run time, and restates the caveat under the decomposition when
you [`print()`](https://rdrr.io/r/base/print.html) the result — that
scan cannot verify the causal structure, and its silence does not
establish that no such confounder exists. VanderWeele & Tchetgen
Tchetgen (2017) propose the randomized interventional analogues
(`mediation_type = "I"`) for that setting.

#### Targeted estimation (TMLE)

For natural effects only, the plug-in simulation can be replaced by the
targeted maximum likelihood estimator of Zheng & van der Laan (2017,
§4.3):

``` r

fit_tmle <- mediation(
  ...,                     # same arguments as above
  mediation_type = "N",
  estimator      = "tmle"  # multiply robust; EIC Wald CIs, no bootstrap
)
```

It runs backward targeted regressions instead of forward simulation.
Zheng & van der Laan (2017) establish multiple robustness for this
estimator — consistency when only certain subsets of the nuisance models
are correct — and it reports Wald CIs from the efficient influence
curve, so `R` and `mc_sample` are ignored. Note that those results are
stated for correctly specified nuisance models, and this implementation
builds its targeted sequential regressions as additive main-effects
working models; see
[`?mediation`](https://adayim.github.io/causalMed/reference/mediation.md)
for the detail. It requires an exposure model, a binary outcome, and
lag-style recodes only (`recodes(lag_A = A)`; exposure lags must be
first-order); derived recodes such as splines or cumulative counts need
`estimator = "gcomp"`. Inspect any positivity warnings it reports before
trusting the affected quantities.

### Enabling parallel bootstrap

``` r

library(future)
plan(multisession)   # Windows; use plan(multicore) on Unix/macOS

fit_parallel <- mediation(..., R = 500)

plan(sequential)     # reset after use
```

## Data format

Input data must be in **long format** with one row per subject per time
point. The time variable should be an ordered integer starting at 0. Lag
variables (e.g. `lag1_A`) must be created by the `init_recode` /
`in_recode` hooks: the Monte Carlo cohort is built from `id_var` and
`base_vars` only, so a lag column that merely exists in the input data
is not carried into the simulation and the run stops with
`object 'lag1_A' not found`. The `nonsurvivaldata` bundled with the
package illustrates the required structure.

``` r

data("nonsurvivaldata", package = "causalMed")
head(nonsurvivaldata)
#>   id time          V         L1 L2 A         M    Y_cont Y_bin lag1_A   lag1_L1
#> 1  1    0  0.4365731  0.7601565  0 1 0.3825810        NA    NA     NA        NA
#> 2  1    1  0.4365731  0.1041467  1 1 1.2635708        NA    NA      1 0.7601565
#> 3  1    2  0.4365731  0.8956876  0 1 0.6277357        NA    NA      1 0.1041467
#> 4  1    3  0.4365731  1.6316564  0 1 1.2583611        NA    NA      1 0.8956876
#> 5  1    4  0.4365731  1.1148361  0 1 0.4602865 0.1888076     1      1 1.6316564
#> 6  2    0 -1.8666578 -0.3374623  0 1 0.3350378        NA    NA     NA        NA
#>   lag1_L2    lag1_M
#> 1      NA        NA
#> 2       0 0.3825810
#> 3       1 1.2635708
#> 4       0 0.6277357
#> 5       0 1.2583611
#> 6      NA        NA
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

7.  Avin, C., Shpitser, I., & Pearl, J. (2005). Identifiability of
    path-specific effects. *Proceedings of the 19th International Joint
    Conference on Artificial Intelligence (IJCAI)*, 357–363.

8.  Miles, C. H. (2023). On the causal interpretation of randomised
    interventional indirect effects. *Journal of the Royal Statistical
    Society: Series B*, 85(4), 1154–1172.
    [doi:10.1093/jrsssb/qkad066](https://doi.org/10.1093/jrsssb/qkad066)
