# Calculate Risk Difference, Risk Ratio, or Mediation Effects based on g-formula

These internal functions compute the Risk Difference (RD) and Risk Ratio
(RR) for g-formula estimates of the total effect, natural direct effect
(NDE), and natural indirect effect (NIE). The functions also calculate
confidence intervals (for RD and RR) and mediation effects based on
g-formula methodology.

## Usage

``` r
risk_estimate_point(data_list, ref_int, intervention, return_data)
```

## Arguments

- data_list:

  A list of simulated datasets returned by `.run_interventions`.

- ref_int:

  Reference intervention to compare against.

- intervention:

  A vector specifying the interventions to compare.

- return_data:

  Logical. If `TRUE`, returns the simulated data along with effects.

## Value

A `data.table` with calculated effects (Risk Difference, Risk Ratio, or
Mediation effects) and, if applicable, confidence intervals. For
mediation effects, the result will show total, direct, and indirect
effects.

## Details

These functions are used internally within the g-formula framework for
effect estimation. They are not intended for direct user input. The
functions perform: - `risk_estimate_point`: Computes point estimates for
risk difference and risk ratio. - `risk_estimate_boot`: Computes
per-bootstrap-replicate estimates for the effects (RD/RR). -
`risk_estimate_mediation`: Computes total, direct, and indirect effects
by subtraction.

These functions operate on the output of the `.run_interventions`
function, estimating effects (both point estimates and confidence
intervals) from simulated data.

## See also

[`gformula`](https://adayim.github.io/causalMed/reference/gformula.md),
[`mediation`](https://adayim.github.io/causalMed/reference/mediation.md)
