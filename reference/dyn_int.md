# Define a dynamic intervention rule

Captures an R expression for a dynamic (rule-based) exposure
intervention. The expression is evaluated at each time step inside the
simulated `data.table`, so any column in the current Monte Carlo dataset
can be referenced by name (including the exposure itself, which holds
its natural-course draw at the time of evaluation).

## Usage

``` r
dyn_int(expr)
```

## Arguments

- expr:

  An R expression that returns a numeric (or logical coerced to numeric)
  vector with one value per row. Column names from the current simulated
  dataset are in scope (e.g. `L1`, `A`).

## Value

An object of class `"causalMed_dynint"`.

## Examples

``` r
# Treat only if the natural-course value of A exceeds 0:
dyn_int(as.numeric(A > 0))
#> [[1]]
#> as.numeric(A > 0)
#> 
#> attr(,"class")
#> [1] "causalMed_dynint"

# Treat if L1 > 0 or L2 equals 1:
dyn_int(as.numeric(L1 > 0 | L2 == 1))
#> [[1]]
#> as.numeric(L1 > 0 | L2 == 1)
#> 
#> attr(,"class")
#> [1] "causalMed_dynint"
```
