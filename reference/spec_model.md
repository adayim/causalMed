# Specify a model for one time-varying variable

Describes the model for one time-varying variable: its formula, its role
in the data-generating process (`mod_type`) and how its values are
simulated (`var_type`). Nothing is fitted here;
[`gformula`](https://adayim.github.io/causalMed/reference/gformula.md)
and
[`mediation`](https://adayim.github.io/causalMed/reference/mediation.md)
fit the models, in list order, and use them to simulate counterfactual
trajectories.

## Usage

``` r
spec_model(
  formula,
  subset = NULL,
  recode = NULL,
  var_type = c("normal", "binary", "categorical", "custom"),
  mod_type = c("covariate", "exposure", "mediator", "outcome", "censor", "survival"),
  custom_fit = NULL,
  custom_sim = NULL,
  truncate = TRUE,
  ...
)
```

## Arguments

- formula:

  Model formula, e.g. `L ~ A + lag1_L + time`. Every variable must exist
  in the analysis data, or be created by a recode hook.

- subset:

  Optional unquoted logical expression, e.g. `platnormm1 == 0`. The
  model is fitted on the rows where it holds and, at each simulated time
  step, draws only those rows; other rows keep their current value. This
  is how an absorbing state is modelled.

- recode:

  Optional
  [`recodes`](https://adayim.github.io/causalMed/reference/recodes.md)
  applied before this model is fitted and before its response is
  simulated.

- var_type:

  How values are simulated:

  `"normal"`

  :   (default) Gaussian draws around the fitted linear model, with the
      residual standard deviation.

  `"binary"`

  :   Bernoulli draws from a logistic model.

  `"categorical"`

  :   Draws from a multinomial logistic model
      ([`multinom`](https://rdrr.io/pkg/nnet/man/multinom.html)),
      returned in the variable's type in the data (a factor keeps its
      levels, numeric codes stay numeric); the Hmisc package must be
      installed.

  `"custom"`

  :   A user-supplied fitting function (`custom_fit`) and/or simulation
      function (`custom_sim`).

  `"censor"` and `"survival"` models must be `"binary"`; an `"outcome"`
  model must be `"binary"` or `"normal"`.

- mod_type:

  Role of the variable: `"covariate"` (default), `"exposure"`,
  `"mediator"`, `"outcome"` (end-of-follow-up outcome), `"survival"`
  (discrete-time event indicator) or `"censor"` (loss to follow-up).
  Under every intervention the censoring indicator is set to zero, so
  risks are the risks under eliminated loss to follow-up, and with
  `estimator = "gcomp"` a censoring model does not change them. It is
  simulated in the natural course of
  [`gformula`](https://adayim.github.io/causalMed/reference/gformula.md)
  and used by `estimator = "tmle"` in
  [`mediation`](https://adayim.github.io/causalMed/reference/mediation.md).

- custom_fit:

  Fitting function for `var_type = "custom"` (ignored, with a warning,
  otherwise); default [`glm`](https://rdrr.io/r/stats/glm.html). Its
  *name* is recorded and evaluated later, so it must be reachable from
  the global environment or a package namespace: a function defined
  inside another function or
  [`local()`](https://rdrr.io/r/base/eval.html) will not be found. A
  namespace-qualified name (e.g.
  [`truncreg::truncreg`](https://rdrr.io/pkg/truncreg/man/truncreg.html))
  also works on parallel bootstrap workers. Without `custom_sim`, the
  fitted object must have a `terms` component and a
  [`coef`](https://rdrr.io/r/stats/coef.html) method.

- custom_sim:

  Simulation function `function(fit, newdata)` returning one simulated
  value per row of `newdata`. It replaces the draw implied by
  `var_type`, so it should return a draw from the model's distribution
  rather than a fitted mean. If it is omitted with
  `var_type = "custom"`, values are drawn from a normal distribution
  around the linear predictor; that is the fitted mean only under an
  identity link, so a fit with another link is rejected. For `"outcome"`
  and `"survival"` models the reported risk is always computed from the
  fitted model, not from `custom_sim`.

- truncate:

  Logical (default `TRUE`). Clip simulated numeric values (`"normal"`
  draws and numeric `custom_sim` output) to the range of the response
  observed in the data. This matches the default `sim_trunc = TRUE` of
  gfoRmula. `FALSE` draws from the untruncated distribution. No effect
  on `"binary"` or `"categorical"` variables.

- ...:

  Further arguments to the fitting function (`glm`, `multinom` or
  `custom_fit`), e.g. `family`.

## Value

An object of class `"causalMed_gmodel"`: the unevaluated fitting call
plus `subset`, `recode`, `var_type`, `mod_type`, `custom_sim` and
`truncate`.

## Details

Choosing a model and a simulation rule that suit each variable is the
analyst's responsibility; the package checks only what it can verify
mechanically.

## See also

[`gformula`](https://adayim.github.io/causalMed/reference/gformula.md),
[`mediation`](https://adayim.github.io/causalMed/reference/mediation.md),
[`recodes`](https://adayim.github.io/causalMed/reference/recodes.md)

## Examples

``` r
# A binary covariate modelled only among those not yet in the state
spec_model(platnorm ~ all + cmv + male + age + gvhdm1 + daysgvhd + wait,
           var_type = "binary", mod_type = "covariate",
           subset = platnormm1 == 0)
#> $call
#> stats::glm(formula = platnorm ~ all + cmv + male + age + gvhdm1 + 
#>     daysgvhd + wait, subset = platnormm1 == 0, family = binomial())
#> 
#> $subset
#> platnormm1 == 0
#> 
#> $recode
#> NULL
#> 
#> $var_type
#> [1] "binary"
#> 
#> $mod_type
#> [1] "covariate"
#> 
#> $custom_sim
#> NULL
#> 
#> $truncate
#> [1] TRUE
#> 
#> attr(,"class")
#> [1] "causalMed_gmodel" "list"            

# A count covariate: Poisson fit with a matching simulation function
sim_poisson <- function(fit, newdata) {
  rpois(nrow(newdata), predict(fit, newdata = newdata, type = "response"))
}
spec_model(daysgvhd ~ all + cmv + male + age + wait,
           var_type = "custom", mod_type = "covariate",
           custom_sim = sim_poisson, family = poisson(link = "log"))
#> $call
#> stats::glm(formula = daysgvhd ~ all + cmv + male + age + wait, 
#>     family = poisson(link = "log"))
#> 
#> $subset
#> NULL
#> 
#> $recode
#> NULL
#> 
#> $var_type
#> [1] "custom"
#> 
#> $mod_type
#> [1] "covariate"
#> 
#> $custom_sim
#> function (fit, newdata) 
#> {
#>     rpois(nrow(newdata), predict(fit, newdata = newdata, type = "response"))
#> }
#> <environment: 0x55a180016690>
#> 
#> $truncate
#> [1] TRUE
#> 
#> attr(,"class")
#> [1] "causalMed_gmodel" "list"            
```
