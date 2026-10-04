# Recode variables inside the simulation

Captures named expressions, unevaluated, to be evaluated inside the
simulated data. They are used for the `init_recode`, `in_recode` and
`out_recode` hooks of
[`gformula`](https://adayim.github.io/causalMed/reference/gformula.md)
and
[`mediation`](https://adayim.github.io/causalMed/reference/mediation.md),
and for the `recode` argument of
[`spec_model`](https://adayim.github.io/causalMed/reference/spec_model.md).
Expressions are evaluated in order, so a later one can use the result of
an earlier one.

The simulated cohort starts with only the subject identifier and the
baseline covariates, so every lag, counter or other derived variable a
model uses must be created here, even if it already exists in the data:

- `init_recode`: at the first time step, before the models are evaluated
  (initial values).

- `in_recode`: at the start of every later time step, before the models
  are evaluated (lags, functions of time).

- `out_recode`: at the end of every time step, the first included, after
  the models are evaluated (cumulative counts, carrying an absorbing
  state forward). A variable it updates needs a starting value from
  `init_recode`.

## Usage

``` r
recodes(...)
```

## Arguments

- ...:

  Named expressions, e.g. `lag1_A = A` or `daysq = day^2`.

## Value

A list of unevaluated expressions of class `"causalMed_recodes"`.

## Examples

``` r
# One-period lags: 0 at the first time step, then the previous value
init <- recodes(lag1_A = 0, lag1_L = 0)
upd  <- recodes(lag1_A = A, lag1_L = L)

# A function of the time variable, e.g. for in_recode
recodes(day_sq = day^2)
#> $day_sq
#> day^2
#> 
#> attr(,"class")
#> [1] "causalMed_recodes" "list"             

# A running count of exposed steps: start at 0 (init_recode), add the
# current step's exposure at the end of each step (out_recode)
init <- recodes(cum_A = 0)
upd  <- recodes(cum_A = cum_A + A)
```
