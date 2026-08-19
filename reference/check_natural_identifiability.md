# Warn about non-identifiability of natural effects under intermediate confounding

For `mediation_type = "N"`, natural direct and indirect effects are not
identifiable from observational data when a confounder of the
mediator-outcome relationship is itself affected by exposure (Avin,
Shpitser & Pearl 2005; VanderWeele 2014; VanderWeele & Tchetgen Tchetgen
2017). This check reads the model formulas only: it reports covariate
models that carry the exposure on the right-hand side, i.e. covariates
the user has modelled as exposure-affected. It does not establish that
such a covariate also confounds the mediator-outcome relationship, and
its silence does not establish that no such confounder exists.

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

Invisibly returns a character vector of the covariate (response) names
whose model includes the exposure on the right-hand side; empty if none.
The caller can persist this so downstream methods (e.g.
`print.gformula`) can re-surface the identifiability caveat.
