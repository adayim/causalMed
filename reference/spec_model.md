# Model specification for G-formula

Specify regression models for the time-varying variables to be used
within the g-formula simulation. This function creates unevaluated
models and can further pass to
[`gformula`](https://adayim.github.io/causalMed/reference/gformula.md)
or
[`mediation`](https://adayim.github.io/causalMed/reference/mediation.md).

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

  an object of class formula: An object of class `formula`: symbolic
  model specification to be fitted (e.g., `Y ~ A + L + time`). The
  variables referenced in `formula` must exist in the analysis
  `data.frame` when fitting/evaluating the model during g-formula
  simulation.

- subset:

  Optional. An unquoted logical expression naming columns of the
  analysis data (e.g. `platnormm1 == 0`) that restricts this model to a
  subset of observations. It is passed through when fitting and
  re-evaluated at each time step during simulation, so only the rows
  satisfying it have their response drawn from this model. Rows that
  never satisfy it keep whatever value they already carry.

- recode:

  Optional. One or more recoding statements built with
  [`recodes`](https://adayim.github.io/causalMed/reference/recodes.md)
  (e.g. `recodes(L_lag1 = L)` or `recodes(M_lag1 = 0)`), applied
  **before** fitting the model and **before** simulating its response
  (useful for dynamic recoding). Anything that is not a
  [`recodes()`](https://adayim.github.io/causalMed/reference/recodes.md)
  object is rejected.

- var_type:

  Character. The response type for simulation/prediction: `"normal"`
  (the default), `"binary"`, `"categorical"`, or `"custom"`. By default,
  values are simulated via:

  - `"binary"`: Bernoulli draws using the fitted mean.

  - `"normal"`: Gaussian draws using fitted mean, **clipped to the
    observed range** of the response (see `truncate`).

  - `"categorical"`: Multinomial draws via
    [`multinom`](https://rdrr.io/pkg/nnet/man/multinom.html).

  - `"custom"`: user-specified via `custom_fit` and/or `custom_sim`
    (numeric output is also clipped to the observed range unless
    `truncate = FALSE`).

- mod_type:

  Character. The role of this model in the data-generating process:
  `"covariate"` (the default), `"exposure"`, `"mediator"`, `"outcome"`,
  `"censor"`, or `"survival"`.

- custom_fit:

  Optional. A model fitting function, used **only** when
  `var_type = "custom"` (it is ignored, with a warning, for the other
  types). If `var_type = "custom"` and `custom_fit` is not provided,
  [`glm`](https://rdrr.io/r/stats/glm.html) is used by default. This can
  be used to define a fitting function other than
  [`glm`](https://rdrr.io/r/stats/glm.html) and
  [`multinom`](https://rdrr.io/pkg/nnet/man/multinom.html).

  **Scope:** `spec_model()` records the function *name*, not the
  function object, and the fit is evaluated inside the package. The name
  must therefore resolve from the global environment or from a package
  namespace — a fitter defined inside another function or inside
  [`local()`](https://rdrr.io/r/base/eval.html) will not be found.
  Prefer a fully qualified name (e.g.
  [`truncreg::truncreg`](https://rdrr.io/pkg/truncreg/man/truncreg.html)
  for truncated regression), which is also what makes the bootstrap work
  under a parallel
  [`plan`](https://future.futureverse.org/reference/plan.html), where
  each worker is a fresh session.

  **What the fitted object must provide:** without `custom_sim` the
  simulation evaluates the fit's linear predictor, so it needs a `terms`
  component and a [`coef`](https://rdrr.io/r/stats/coef.html) method; a
  fit lacking either is rejected with an explanatory error naming the
  variable. With `custom_sim` the drawing is delegated, and neither is
  required. Choosing a model that is appropriate for the variable, and a
  `custom_sim` that draws from the distribution that model implies,
  remains the analyst's responsibility.

- custom_sim:

  Optional. A simulation function for the model. It must accept two
  arguments – the fitted model object and a `data.frame` of new data to
  predict on – and return a vector of simulated responses of matching
  length. When supplied it takes priority over the drawing rule implied
  by `var_type`. If omitted and `var_type = "custom"`, normal draws are
  used by default, which requires the fitted object to have a `terms`
  component and a [`coef`](https://rdrr.io/r/stats/coef.html) method.

  **It supplies a draw, not a fitted value.** At each time step the
  Monte Carlo engine assigns this variable the vector `custom_sim`
  returns, in place of the draw the built-in `var_type` rule would have
  made (a Bernoulli draw for `"binary"`, a Gaussian draw around the
  linear predictor for `"normal"`). Returning
  [`predict()`](https://rdrr.io/r/stats/predict.html)'s fitted value
  therefore assigns the conditional mean rather than a draw from the
  conditional distribution. Whether that is appropriate for the variable
  being simulated is the analyst's decision.

  **It does not apply to the outcome.** For `mod_type = "outcome"` or
  `"survival"` the reported risk is computed from the fitted model's own
  linear predictor, so `custom_sim` affects only the simulated response
  value, never the estimate. Those models must therefore have
  extractable coefficients.

- truncate:

  Logical. If `TRUE` (default), simulated numeric values are clipped to
  the range of the response observed in the data — i.e. a `"normal"`
  draw is clipped to `[min, max]` of the observed response, and numeric
  output of `custom_sim` is clipped as well. `TRUE` is the default. It
  corresponds to the `sim_trunc` argument of gfoRmula, documented there
  as "whether to truncate simulated covariates to their range in the
  observed data set", whose default is also `TRUE`. Set `FALSE` to draw
  from the untruncated fitted distribution (corresponding to
  `sim_trunc = FALSE`), or to let a `custom_sim` function be
  authoritative over its own output range. Which is appropriate depends
  on the variable being simulated and is the analyst's decision. Has no
  effect on `"binary"` or `"categorical"` responses.

- ...:

  Other parameters passed to the model fitting function,
  [`glm`](https://rdrr.io/r/stats/glm.html),
  [`multinom`](https://rdrr.io/pkg/nnet/man/multinom.html) or
  `custom_fit`.

## Value

An object of class `"causalMed_gmodel"`:

- `call`:

  An unevaluated call (function + arguments) to fit the model.

- `subset`:

  The unevaluated subset expression provided via `subset`.

- `recode`:

  The recoding statements provided via `recode`.

- `var_type`:

  The response type, as provided.

- `mod_type`:

  The model role, as provided.

- `custom_sim`:

  The simulation function provided via `custom_sim`.

- `truncate`:

  The truncation flag, as provided.

## Details

This function will be used to create an unevaluated model for the
g-formula. `spec_model()` does not fit the model immediately. It returns
an unevaluated call plus metadata (`var_type`, `mod_type`, `subset`,
`recode`, `custom_sim`) that are used later by
[`gformula`](https://adayim.github.io/causalMed/reference/gformula.md)/[`mediation`](https://adayim.github.io/causalMed/reference/mediation.md)
to fit in temporal order and to simulate counterfactual trajectories.

## See also

[`gformula`](https://adayim.github.io/causalMed/reference/gformula.md),[`mediation`](https://adayim.github.io/causalMed/reference/mediation.md)

## Examples

``` r
data(gvhd)
mod_cov1 <- spec_model(platnorm ~ all + cmv + male + age + agecurs1 +
  agecurs2 + gvhdm1 + daysgvhd + daysnorelapse + wait,
var_type = "binary",
mod_type = "covariate",
subset = platnormm1 == 0
)
## For Poisson regression
predict_poisson <- function(fit, newdf) {
  theta <- stats::predict(object = fit, type = "response", newdata = newdf)
  prediction <- rpois(n = nrow(newdf), lambda = theta)
  return(prediction)
}
mod_cov1 <- spec_model(platnorm ~ all + cmv + male + age + agecurs1 +
  agecurs2 + gvhdm1 + daysgvhd + daysnorelapse + wait,
var_type = "custom",
mod_type = "covariate",
subset = platnormm1 == 0,
custom_sim = predict_poisson,
family = "poisson"(link = "log"),
y = TRUE
)
```
