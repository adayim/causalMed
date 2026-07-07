# Example Dataset for a Non-Survival Outcome

A simulated dataset with time-varying and baseline variables for 1000
subjects over 5 time points, including exposure, mediator, confounders,
and outcome.

## Usage

``` r
nonsurvivaldata
```

## Format

A data frame with 5000 rows and 13 variables:

- id:

  Unique subject identifier.

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

  Time-varying mediator.

- Y_bin:

  Binary outcome observed at each time point.

- Y_cont:

  Continuous outcome observed at each time point.

- lag1_A:

  Lagged exposure (A at previous time point).

- lag1_L1:

  Lagged confounder L1.

- lag1_L2:

  Lagged confounder L2.

- lag1_M:

  Lagged mediator.

## Source

Simulated data generated for package examples.

## Details

The simulated longitudinal data-generating structure can be summarized
as: \$\$ A_t \leftarrow V, L1\_{t-1}, L2\_{t-1}, A\_{t-1}, t;\quad L1_t
\leftarrow V, A_t, L1\_{t-1}, t;\quad L2_t \leftarrow V, A_t, L2\_{t-1},
t;\quad M_t \leftarrow V, A_t, L1_t, L2_t, M\_{t-1}, t;\quad Y_t
\leftarrow V, A_t, M_t, L1_t, L2_t, A_t\*M_t. \$\$ The same outcome
model structure is used for both `Y_bin` and `Y_cont`, with the
appropriate outcome distribution specified for each outcome type.
