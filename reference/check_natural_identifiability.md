# Warn about non-identifiability of natural effects under intermediate confounding

For `mediation_type = "N"`, natural direct and indirect effects are not
identifiable from observational data when an intermediate confounder of
the mediator-outcome relationship is itself affected by exposure (Avin,
Shpitser & Pearl 2005; VanderWeele 2014; VanderWeele & Tchetgen Tchetgen
2017). This check inspects the covariate models for dependence on the
exposure variable on the right-hand side of the formula; if any are
found, a warning is emitted suggesting `mediation_type = "I"` instead.

## Usage

``` r
check_natural_identifiability(models, exposure)
```

## Arguments

- models:

  List of model specifications from
  [`spec_model`](https://adayim.github.io/causalMed/reference/spec_model.md).

- exposure:

  Character scalar. Name of the exposure variable.

## Value

Invisibly returns `NULL`.
