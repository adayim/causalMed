# Example Dataset for a Survival Outcome

A simulated dataset with time-varying and baseline variables for 3000
subjects with a discrete-time survival outcome, suitable for use with
the g-formula and mediation functions under survival settings.

## Usage

``` r
survivaldata
```

## Format

A data frame with 7113 rows and 10 variables, in long at-risk format
(one row per subject per time point while still at risk; 3000 subjects,
2799 events, no missing values):

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

  Time-varying mediator (continuous).

- Y:

  Event indicator at this time point (1 = event, 0 = still at risk).
  There is no censoring in this dataset.

- lag1_A:

  Exposure at the previous time point; 0 at `time == 0`.

- lag1_M:

  Previous value of `M`; 0 at `time == 0`.

- lag1_L:

  Previous value of `L`; 0 at `time == 0`.

## Source

Simulated data generated for package examples.

## Details

The data-generating structure can be summarized as: \$\$ A_t \leftarrow
V, L\_{t-1}, A\_{t-1}, t;\quad L_t \leftarrow V, A_t, L\_{t-1}, t;\quad
M_t \leftarrow V, A_t, L_t, M\_{t-1}, t;\quad Y_t \leftarrow V, A_t,
M_t, L_t, A_t\*M_t, t. \$\$

Events and follow-up observations are generated only while subjects
remain at risk. Therefore, once `Y` becomes 1 at time \\t\\, no
observations are retained for that subject at subsequent time points
\\t+1, t+2, \ldots\\.

## See also

[`nonsurvivaldata`](https://adayim.github.io/causalMed/reference/nonsurvivaldata.md)
for the end-of-follow-up counterpart.
