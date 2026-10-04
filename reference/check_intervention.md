# Check for the intervention

Check if the intervention is correctly defined.

## Usage

``` r
check_intervention(models, intervention, ref_int, time_len)
```

## Arguments

- models:

  List of
  [`spec_model`](https://adayim.github.io/causalMed/reference/spec_model.md)
  objects, in the order the variables are generated.

- intervention:

  Named list of interventions. Each element is `NULL` (the natural
  course: exposure drawn from its fitted model), a 0/1 value or vector
  with one value per distinct time point (a static intervention), or a
  [`dyn_int`](https://adayim.github.io/causalMed/reference/dyn_int.md)
  rule, e.g.
  `list(natural = NULL, treat_if_high = dyn_int(as.numeric(L1 > 0)))`.
  `NULL` (the default) runs the natural course only.

- ref_int:

  Reference for the contrasts: `0` or `"natural"` (default) for the
  natural course, or the position or name of an element of
  `intervention`. See Details.

- time_len:

  length of the time in the data.
