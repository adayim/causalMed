# Example Dataset for a Non-Survival Outcome

A simulated dataset with time-varying and baseline variables for 3000
subjects over 5 time points, including exposure, mediator, confounders,
and an end-of-follow-up outcome. Every subject contributes all 5 rows
(there is no attrition), which makes this the simpler of the two
simulated examples.

## Usage

``` r
nonsurvivaldata
```

## Format

A data frame with 15,000 rows (3000 subjects x 5 time points) and 13
variables:

- id:

  Unique subject identifier (1-3000).

- time:

  Time variable (0 to 4).

- V:

  Time-fixed baseline covariate.

- L1:

  Time-varying confounder 1 (continuous).

- L2:

  Time-varying confounder 2 (binary).

- A:

  Time-varying binary exposure.

- M:

  Time-varying mediator (continuous).

- Y_bin:

  Binary end-of-follow-up outcome. Recorded **only at `time == 4`**;
  `NA` at `time` 0-3.

- Y_cont:

  Continuous end-of-follow-up outcome, on the same schedule as `Y_bin`:
  recorded **only at `time == 4`**, `NA` otherwise.

- lag1_A:

  Exposure at the previous time point; `NA` at `time == 0`.

- lag1_L1:

  Previous value of `L1`; `NA` at `time == 0`.

- lag1_L2:

  Previous value of `L2`; `NA` at `time == 0`.

- lag1_M:

  Previous value of `M`; `NA` at `time == 0`.

## Source

Simulated from the parametric process given under Details.

## Details

The simulated longitudinal data-generating structure can be summarized
as: \$\$ A_t \leftarrow V, L1\_{t-1}, L2\_{t-1}, A\_{t-1}, t;\quad L1_t
\leftarrow V, A_t, L1\_{t-1}, t;\quad L2_t \leftarrow V, A_t, L2\_{t-1},
t;\quad M_t \leftarrow V, A_t, L1_t, L2_t, M\_{t-1}, t;\quad Y
\leftarrow V, A_4, M_4, L1_4, L2_4, A_4\*M_4. \$\$ The exposure,
mediator and confounders evolve at every time point, but the outcome is
realised once, at the end of follow-up.

The generating parameters are, with \\V \sim N(0, 1)\\ and all lagged
terms set to zero at \\t = 0\\:


      logit P(A_t = 1) =  1.25 + 0.50 V + 0.25 L1_{t-1} + 0.25 L2_{t-1}
                               + 0.20 A_{t-1} + 0.01 t
      E(L1_t)          =  0.02 + 0.50 V + 0.18 A_t + 0.35 L1_{t-1} + 0.02 t
                         (Gaussian, SD 0.70)
      logit P(L2_t = 1)= -0.10 + 0.40 V + 0.18 A_t + 0.35 L2_{t-1} + 0.02 t
      E(M_t)           =  0.05 + 0.50 V + 0.45 A_t + 0.12 L1_t + 0.10 L2_t
                               + 0.30 M_{t-1} + 0.02 t   (Gaussian, SD 0.50)
      logit(pi)        = -2.60 + 0.50 V + 0.51 A_4 + 0.45 M_4
                               + 0.25 L1_4 - 0.25 L2_4 + 0.30 A_4 M_4

Both outcomes share that same \\\pi\\: `Y_bin` is
\\\mathrm{Bernoulli}(\pi)\\ and `Y_cont` is \\\mathrm{Beta}(50\pi,\\
50(1-\pi))\\. `Y_cont` therefore lies strictly in \\(0, 1)\\; fitting it
with `var_type = "normal"` is a working approximation, not the
generating distribution.

The confounders \\L1\\ and \\L2\\ depend on current exposure, so
exposure-affected confounding of the mediator-outcome relationship is
present – the setting the interventional effects are designed for.

**Large-sample true values** (risk differences, computed under the same
permuted-pool construction the package uses for `mediation_type = "I"`,
at 200,000 subjects):

|                                           |       |
|-------------------------------------------|-------|
| Interventional direct effect (IDE)        | 0.075 |
| Interventional indirect effect (IIE)      | 0.071 |
| Total effect (TE)                         | 0.155 |
| Decomposition residual (TE \\-\\ overall) | 0.009 |

The two pathways are deliberately of comparable size, and the residual
is deliberately non-zero, so that the decomposition and its
non-additivity are both visible in a worked example.

Two consequences for model specification:

- Because the outcome exists at a single time point, `time` is constant
  among the rows used to fit the outcome model and must **not** appear
  in its formula: the term is not estimable and would be dropped (with a
  warning) from the simulation.

- The `lag1_*` columns are `NA` at `time == 0`, so a baseline value must
  be supplied through `init_recode`, e.g.
  `init_recode = recodes(lag1_A = 0, lag1_L1 = 0)`. (The lag columns of
  [`survivaldata`](https://adayim.github.io/causalMed/reference/survivaldata.md)
  are 0-filled at baseline instead, so no initialisation is strictly
  required there.)

## See also

[`survivaldata`](https://adayim.github.io/causalMed/reference/survivaldata.md)
for the survival-outcome counterpart.
