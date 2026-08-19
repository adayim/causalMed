# Check that a natural-course intervention can be simulated

A natural-course intervention – a `NULL` element of `intervention`, or
`intervention = NULL` – draws the exposure from its own fitted model at
each time step. Without a `mod_type = "exposure"` model there is nothing
to draw from, and the exposure column never enters the Monte Carlo
dataset, so the run fails deep inside
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) with an
opaque "object '\<exposure\>' not found". Fail early with an actionable
message.

## Usage

``` r
check_natural_course(models, intervention, exposure)
```

## Arguments

- models:

  List of model specifications from
  [`spec_model`](https://adayim.github.io/causalMed/reference/spec_model.md).

- intervention:

  The resolved (named, non-`NULL`) intervention list.

- exposure:

  Character scalar. Name of the exposure variable.
