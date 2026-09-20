# Validate one static intervention vector

One static exposure regime is a numeric or logical vector, with no `NA`,
of length 1 (the same value at every time point) or `time_len` (one
value per **distinct** time point in the data, in sorted order). Shared
by
[`check_intervention`](https://adayim.github.io/causalMed/reference/check_intervention.md)
(for
[`gformula()`](https://adayim.github.io/causalMed/reference/gformula.md))
and
[`mediation`](https://adayim.github.io/causalMed/reference/mediation.md)
(for its two exposure regimes).

## Usage

``` r
check_static_regime(x, time_len, arg, values = NULL)
```

## Arguments

- x:

  The vector to check.

- time_len:

  Number of distinct time points in the data.

- arg:

  Argument name used in error messages.

- values:

  Optional vector of allowed values. `NULL` (the default) allows any
  numeric value.

## Value

`x` coerced to numeric.
