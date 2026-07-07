# Check temporal ordering of models for mediation analysis

Issues a warning if covariate or exposure models appear in positions
that violate the assumed A(t) -\> M(t) -\> L(t) -\> S(t) ordering.

## Usage

``` r
check_mediation_order(models)
```

## Arguments

- models:

  List of model specifications from
  [`spec_model`](https://adayim.github.io/causalMed/reference/spec_model.md).
