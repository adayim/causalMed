# Example Dataset for a Survival Outcome

A simulated dataset with time-varying and baseline variables for
subjects with a survival outcome, suitable for use with the g-formula
and mediation functions under survival settings.

## Usage

``` r
survivaldata
```

## Format

A data frame with 7113 rows and 10 variables:

- id:

  Unique subject identifier.

- time:

  Time variable (integer, starting at 0).

- V:

  Time-fixed baseline covariate.

- L:

  Time-varying confounder.

- A:

  Time-varying binary exposure.

- M:

  Time-varying mediator.

- Y:

  Survival outcome indicator (1 = event, 0 = alive/censored).

- lag1_A:

  Lagged exposure (A at previous time point).

- lag1_M:

  Lagged mediator.

- lag1_L:

  Lagged confounder.

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
