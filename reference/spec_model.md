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

  Optional. A logical expression or character vector defining a subset
  of observations to fit/simulate this model on; passed through when
  fitting and also respected during simulation of the response for this
  model.

- recode:

  Optional character vector. Should be defined with
  [`recodes`](https://adayim.github.io/causalMed/reference/recodes.md).
  One or more recoding statements (e.g., `L_lag1 = L` or `M_lag1 = 0`)
  to be applied \*\*before\*\* evaluating the model and \*\*before\*\*
  simulating the response for this model (useful for dynamic recoding).

- var_type:

  Character. The response type for simulation/prediction: `"binary"`,
  `"normal"`, `"categorical"`, or `"custom"`. By default, values are
  simulated via:

  - `"binary"`: Bernoulli draws using the fitted mean.

  - `"normal"`: Gaussian draws using fitted mean.

  - `"categorical"`: Multinomial draws via
    [`multinom`](https://rdrr.io/pkg/nnet/man/multinom.html).

  - `"custom"`: user-specified via `custom_fit` and/or `custom_sim`.

- mod_type:

  Character. The role of this model in the data-generating process:
  `"covariate"`, `"exposure"`, `"mediator"`, `"outcome"`, `"censor"`, or
  `"survival"`.

- custom_fit:

  Optional. A model fitting function used when `var_type = "custom"`. If
  `var_type = "custom"` and `custom_fit` is not provided,
  [`glm`](https://rdrr.io/r/stats/glm.html) is used by default. This can
  be used to define the modeling fitting function other than
  [`glm`](https://rdrr.io/r/stats/glm.html) and
  [`multinom`](https://rdrr.io/pkg/nnet/man/multinom.html). For parallel
  computation, provide the fully qualified function name (e.g.,
  [`truncreg::truncreg`](https://rdrr.io/pkg/truncreg/man/truncreg.html)
  for truncated regression).

- custom_sim:

  Optional. A simulation function for the model. It should accept
  exactly two arguments: the fitted model object and a `data.frame` of
  new data to predict on, and return a vector of simulated responses of
  matching length. If omitted and `var_type = "custom"`, normal draws
  are used by default.

- ...:

  Other parameters passed to the model fitting function,
  [`glm`](https://rdrr.io/r/stats/glm.html),
  [`multinom`](https://rdrr.io/pkg/nnet/man/multinom.html) or
  `custom_fit`.

## Value

An object of class `"gmodel"`:

- `call`:

  An unevaluated call (function + arguments) to fit the model.

- `subset`:

  The subset expression/character vector provided via `subset`.

- `recode`:

  The recoding statements provided via `recode`.

- `var_type`:

  The response type, as provided.

- `mod_type`:

  The model role, as provided.

- `custom_sim`:

  A symbolic customized simulation function name, as provided.

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
  theta <- stats::predict(object = fitcov, type = "response", newdata = newdf)
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
