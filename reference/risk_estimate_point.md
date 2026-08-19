# Point-estimate contrasts against a reference intervention

Internal use only. Turns the per-intervention mean outcomes produced by
`.run_interventions` into risk differences and risk ratios against the
reference intervention. Confidence intervals are *not* computed here;
[`gformula`](https://adayim.github.io/causalMed/reference/gformula.md)
adds them from the bootstrap replicates.

## Usage

``` r
risk_estimate_point(data_list, ref_int, intervention, return_data)
```

## Arguments

- data_list:

  Named list of intervention results from `.run_interventions`: one mean
  outcome per intervention, or one simulated `data.table` per
  intervention when `return_data` is `TRUE`.

- ref_int:

  Character scalar naming the reference element of `data_list`.

- intervention:

  The named intervention list; its names determine which contrasts are
  formed.

- return_data:

  Logical. `TRUE` when `data_list` holds simulated datasets rather than
  scalar means.

## Value

A `data.table` with columns `Intervention`, `Risk_type` (`"Difference"`
or `"Ratio"`) and `Estimate`: two rows per non-reference intervention.

## Details

Companion internal functions in the same file, not documented
separately: `risk_estimate_boot` does the same job for one stacked
bootstrap replicate, and `risk_estimate_mediation` builds the mediation
decomposition (total, direct, per-mediator indirect) by subtraction.

## See also

[`gformula`](https://adayim.github.io/causalMed/reference/gformula.md),
[`mediation`](https://adayim.github.io/causalMed/reference/mediation.md)
