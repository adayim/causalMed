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

Simulated data generated for package examples.

## Details

The simulated longitudinal data-generating structure can be summarized
as: \$\$ A_t \leftarrow V, L1\_{t-1}, L2\_{t-1}, A\_{t-1}, t;\quad L1_t
\leftarrow V, A_t, L1\_{t-1}, t;\quad L2_t \leftarrow V, A_t, L2\_{t-1},
t;\quad M_t \leftarrow V, A_t, L1_t, L2_t, M\_{t-1}, t;\quad Y
\leftarrow V, A_4, M_4, L1_4, L2_4, A_4\*M_4. \$\$ The exposure,
mediator and confounders evolve at every time point, but the outcome is
realised once, at the end of follow-up. The same outcome structure is
used for both `Y_bin` and `Y_cont`, with the appropriate distribution
specified for each.

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
