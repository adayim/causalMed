# Check that the natural-effect mediator can be evaluated on the other regime

Under `mediation_type = "N"` the cross-world mediator is drawn from its
model evaluated on the intervention's own covariate history with the
exposure history set to the other regime (Zheng & van der Laan 2017, Eq.
5). The Monte Carlo engine sets exactly two kinds of input to that
regime: the exposure itself (so exposure terms written in the mediator
formula, e.g. `A:L`, are evaluated on it) and first-order exposure lags
(an `in_recode` entry that only copies the exposure, e.g.
`recodes(lag_A = A)`). Every other column a recode derives from the
exposure – a chained lag, a cumulative count or other expression, an
`out_recode` copy, a column created by a model's own `recode` – keeps
the intervention's own exposure history, and the mediator `subset` is
evaluated on the intervention's own data. The check fails closed: a
mediator formula that reads such a column, or a mediator `subset` that
reads the exposure or anything derived from it, is rejected. A mediator
`custom_sim` is handed the whole data set rather than the formula's
variables, so it is rejected too whenever any such column exists.

## Usage

``` r
check_natural_exposure_history(
  models,
  exposure,
  init_recode = NULL,
  in_recode = NULL,
  out_recode = NULL
)
```

## Arguments

- models:

  List of model specifications from
  [`spec_model`](https://adayim.github.io/causalMed/reference/spec_model.md).

- exposure:

  Character scalar. Name of the exposure variable.

- init_recode, in_recode, out_recode:

  The recode hooks passed to
  [`mediation`](https://adayim.github.io/causalMed/reference/mediation.md).

## Value

Invisibly, the names of the first-order exposure lag columns.
