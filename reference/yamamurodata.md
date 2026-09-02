# Yamamuro et al. (2021) Multi-Mediator Simulation Dataset

A single dataset simulated from the data-generating process of the
Yamamuro et al. (2021) multi-mediator simulation study: a time-varying
binary treatment, a time-varying confounder, two sequential continuous
mediators, and a discrete-time survival outcome over three visits. The
large-sample true interventional effects for this process are documented
below, so estimates obtained with
[`mediation`](https://adayim.github.io/causalMed/reference/mediation.md)
can be compared against known values.

## Usage

``` r
yamamurodata
```

## Format

A data frame with 29482 rows and 15 variables (long at-risk format: one
row per subject per visit until the event):

- id:

  Unique subject identifier.

- time:

  Visit index (0, 1, 2).

- V:

  Time-fixed binary baseline covariate.

- A:

  Time-varying binary treatment.

- L:

  Time-varying continuous confounder.

- M1:

  First mediator (continuous), affected by `A` and `L`.

- M2:

  Second mediator (continuous), affected by `A`, `L`, and `M1`.

- Y:

  Event indicator (1 = event at this visit).

- lag1_A, lag1_L, lag1_M1, lag1_M2:

  Previous-visit values (0 at baseline).

- L0base, M10base, M20base:

  Baseline (visit 0) values of `L`, `M1`, `M2`, carried as time-fixed
  columns so the Monte Carlo simulation can be seeded via `init_recode`.

## Source

Simulated from the data-generating process of Yamamuro, S., Shinozaki,
T., Iimuro, S., & Matsuyama, Y. (2021). Mediational g-formula for
time-varying treatment and repeated-measured multiple mediators.
*Statistical Methods in Medical Research*, 30(8), 1782-1799.
[doi:10.1177/09622802211025988](https://doi.org/10.1177/09622802211025988)

## Details

The within-visit ordering is **A \\\to\\ L \\\to\\ M1 \\\to\\ M2 \\\to\\
Y**; the full data-generating equations are translated from the SAS
`%simdata` macro in the paper's supplementary material and are
reproduced in `data-raw/yamamurodata.R` in the package source
repository.

The published study design generates 1000 replicate datasets of \\n =
1000\\ subjects and averages the estimates. A single replicate of that
size is dominated by sampling error (between-replicate SD of the total
effect is about 1.7 percentage points), so the dataset shipped here is
one draw of \\n = 10000\\ subjects (seed 2468) from the identical
process, making single-dataset estimates informative.

**Large-sample true values** (risk differences in percentage points,
computed from the data-generating process at \\n = 10^7\\; Monte Carlo
SE in parentheses):

|                                             |                   |
|---------------------------------------------|-------------------|
| Total effect (TE = \\E\[Y_1\] - E\[Y_0\]\\) | \\-6.36\\ (0.011) |
| Interventional direct effect (IDE)          | \\-3.20\\ (0.013) |
| Interventional indirect effect via M1       | \\-2.29\\ (0.011) |
| Interventional indirect effect via M2       | \\-0.97\\ (0.009) |
| Decomposition residual (TE \\-\\ overall)   | \\0.10\\ (0.016)  |

The decomposition residual is non-zero because the reported total effect
is a natural-course contrast while the direct and indirect effects are
interventional; see the *Mediator pool* section of
[`mediation`](https://adayim.github.io/causalMed/reference/mediation.md)
for how the mediator draws are constructed.
