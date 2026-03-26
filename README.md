
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

# Specify models in temporal order: A -> L -> M -> Y
models_total <- list(
  spec_model(A ~ V + A_lag1 + L_lag1 + time,
             var_type = "binary",  mod_type = "exposure"),
  spec_model(L ~ V + A + L_lag1 + time,
             var_type = "normal",  mod_type = "covariate"),
  spec_model(M ~ V + A + L + M_lag1 + time,
             var_type = "normal",  mod_type = "mediator"),
  spec_model(Y ~ V + A + M + L + A:M,
             var_type = "binary",  mod_type = "outcome")
)

fit_total <- gformula(
  data         = nonsurvivaldata,
  id_var       = "id",
  time_var     = "time",
  base_vars    = "V",
  exposure     = "A",
  models       = models_total,
  intervention = list(never = 0, always = 1),
  ref_int      = "never",
  init_recode  = recodes(A_lag1 = 0, L_lag1 = 0, M_lag1 = 0),
  in_recode    = recodes(A_lag1 = A, L_lag1 = L, M_lag1 = M),
  mc_sample    = 10000,
  R            = 200,
  seed         = 12345
)

print(fit_total)
```

### Mediation analysis: interventional IDE/IIE (Lin et al. 2017)

``` r
library(causalMed)

data("nonsurvivaldata", package = "causalMed")

# Model list must include a mediator model (mod_type = "mediator")
# List order must reflect the temporal DAG: A -> L -> M -> Y
models_med <- list(
  spec_model(A ~ V + A_lag1 + L_lag1 + time,
             var_type = "binary",  mod_type = "exposure"),
  spec_model(L ~ V + A + L_lag1 + time,
             var_type = "normal",  mod_type = "covariate"),
  spec_model(M ~ V + A + L + M_lag1 + time,
             var_type = "normal",  mod_type = "mediator"),
  spec_model(Y ~ V + A + M + L + A:M,
             var_type = "binary",  mod_type = "outcome")
)

fit_med <- mediation(
  data           = nonsurvivaldata,
  id_var         = "id",
  time_var       = "time",
  base_vars      = "V",
  exposure       = "A",
  outcome        = "Y",
  models         = models_med,
  init_recode    = recodes(A_lag1 = 0, L_lag1 = 0, M_lag1 = 0),
  in_recode      = recodes(A_lag1 = A, L_lag1 = L, M_lag1 = M),
  mediation_type = "I",   # interventional IDE/IIE (Lin et al. 2017)
  mc_sample      = 10000,
  R              = 200,   # bootstrap replicates; set R = 1 for point estimate only
  seed           = 12345
)

print(fit_med)
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
