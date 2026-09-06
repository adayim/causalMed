# Simulate Data

Loop through the models and apply any recoding or subset.

## Usage

``` r
simulate_data(
  data,
  exposure,
  models,
  intervention = NULL,
  mediation_type = c(NA, "N", "I"),
  med_pool = NULL
)
```

## Arguments

- data:

  Data to be used for the data generation

- models:

  Model list passed from
  [`gformula`](https://adayim.github.io/causalMed/reference/gformula.md)
  or
  [`mediation`](https://adayim.github.io/causalMed/reference/mediation.md).

- intervention:

  One of:

  - `NULL` — natural-course draw of the exposure (gformula).

  - Numeric scalar/vector or `causalMed_dynint` — static or dynamic
    exposure rule for
    [`gformula()`](https://adayim.github.io/causalMed/reference/gformula.md).

  - `causalMed_intervention` — structured mediation intervention
    carrying a fixed treatment value plus an optional named list of
    per-mediator overrides. See `intervention_spec()`.

- mediation_type:

  Type of the mediation analysis, if the value is `NA` no mediation
  analysis will be performed (default).

- med_pool:

  Optional named list keyed by mediator response variable. Each element
  is the time-\\t\\ slice of the pre-permuted joint mediator-trajectory
  pool built by `.run_interventions` from the corresponding reference
  intervention (Phi00 or Phi11). When the intervention `intervention`
  requires a mediator override under `mediation_type = "I"`, the
  mediator is assigned directly from this vector. This is the joint,
  whole-population marginal draw of the Lin et al. (2017, *Stat Med*)
  Section 4 algorithm and the reference SAS macros (mGFORMULA; Yamamuro
  et al. 2021 Figure 3 step 3). Their Eq. 4 and Eq. 2 are written
  conditional on baseline covariates; see the Mediator pool section of
  [`mediation`](https://adayim.github.io/causalMed/reference/mediation.md).
