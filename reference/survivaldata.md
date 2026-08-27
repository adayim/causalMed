# Example Dataset for a Survival Outcome

A simulated dataset with time-varying and baseline variables for 3000
subjects with a discrete-time survival outcome, suitable for use with
the g-formula and mediation functions under survival settings.

## Usage

``` r
survivaldata
```

## Format

A data frame with 10,223 rows and 11 variables, in long at-risk format
(one row per subject per time point while still at risk; 3000 subjects,
1822 events, 283 censored, no missing values):

- id:

  Unique subject identifier (1-3000).

- time:

  Time variable (integer, 0 to 4).

- V:

  Time-fixed baseline covariate.

- L:

  Time-varying confounder (continuous).

- A:

  Time-varying binary exposure.

- M:

  Time-varying binary mediator.

- Y:

  Event indicator at this time point (1 = event, 0 = still at risk).

- C:

  Loss-to-follow-up indicator at this time point (1 = censored).
  Censoring is drawn only among subjects who did not have the event at
  that time, so `Y` and `C` are never both 1.

- lag1_A:

  Exposure at the previous time point; 0 at `time == 0`.

- lag1_M:

  Previous value of `M`; 0 at `time == 0`.

- lag1_L:

  Previous value of `L`; 0 at `time == 0`.

## Source

Simulated from the parametric process given under Details.

## Details

The data-generating structure can be summarized as: \$\$ A_t \leftarrow
V, L\_{t-1}, A\_{t-1}, t;\quad L_t \leftarrow V, A_t, L\_{t-1}, t;\quad
M_t \leftarrow V, A_t, L_t, M\_{t-1}, t;\quad Y_t \leftarrow V, A_t,
M_t, L_t, A_t\*M_t, t. \$\$

The generating parameters are, with \\V \sim N(0, 1)\\ and all lagged
terms set to zero at \\t = 0\\:


      logit P(A_t = 1) = -0.40 + 0.60 V + 0.40 L_{t-1} + 0.58 A_{t-1} + 0.05 t
      E(L_t)           =  0.50 + 0.60 V + 0.30 A_t + 0.20 L_{t-1} + 0.05 t
                          (Gaussian, SD 0.10)
      logit P(M_t = 1) = -0.95 + 0.50 V + 1.50 A_t + 0.40 L_t
                               + 0.20 M_{t-1} + 0.05 t
      logit P(Y_t = 1) = -3.40 - 0.50 V + 0.25 A_t + 1.70 M_t + 0.30 L_t
                               + 0.50 A_t M_t
      logit P(C_t = 1) = -3.80 + 0.30 V + 0.30 L_t + 0.10 t

The hazard carries no time trend. The confounder \\L\\ depends on
current exposure, so exposure-affected confounding of the
mediator-outcome relationship is present.

**Large-sample true values** (cumulative risk by the fifth period,
computed under the same permuted-pool construction the package uses for
`mediation_type = "I"`, at 200,000 subjects):

|                                                  |       |
|--------------------------------------------------|-------|
| Risk under no exposure                           | 0.39  |
| Risk under exposure                              | 0.79  |
| Interventional direct effect (IDE)               | 0.194 |
| Interventional indirect effect (IIE)             | 0.191 |
| Total effect (TE)                                | 0.400 |
| Mediated-interaction residual (TE \\-\\ overall) | 0.015 |

Follow-up stops at the event or at loss to follow-up, whichever comes
first, so no rows are retained for a subject after either. Of the 3000
subjects, 1822 have the event, 283 are censored and the remaining 895
reach the end of the fifth period.

Censoring depends on \\V\\ and on the current confounder \\L\\, so it is
informative given the observed history and a censoring model is required
for valid estimation: pass one as a `mod_type = "censor"` entry in
`models`. Because censoring is independent of the counterfactual event
process given that history, the target estimand – risk with censoring
eliminated – is unchanged by its presence, so the true values below are
the same as they would be without it.

## See also

[`nonsurvivaldata`](https://adayim.github.io/causalMed/reference/nonsurvivaldata.md)
for the end-of-follow-up counterpart.
