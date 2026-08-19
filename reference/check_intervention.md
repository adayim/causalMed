# Check for the intervention

Check if the intervention is correctly defined.

## Usage

``` r
check_intervention(models, intervention, ref_int, time_len)
```

## Arguments

- models:

  A list of model specifications evaluated in temporal order. The order
  appeared in the list should reflect the temporal ordering of the
  variables, in another way data generation process. See
  [`spec_model`](https://adayim.github.io/causalMed/reference/spec_model.md)
  for a recommended constructor.

- intervention:

  A named list specifying exposure interventions. Each element is one
  of:

  - `NULL` — the natural course (exposure drawn from its fitted model).

  - A numeric/logical scalar or vector (length 1 or equal to the number
    of time points) — a static intervention setting the exposure to that
    value at every (or each specific) time step.

  - A
    [`dyn_int`](https://adayim.github.io/causalMed/reference/dyn_int.md)
    object — a dynamic (rule-based) intervention whose expression is
    evaluated inside the simulated dataset at each time step. Column
    names (including the exposure after its natural-course draw) are in
    scope, e.g.
    `list(natural = NULL, threshold = dyn_int(as.numeric(A > 0)))`.

  If `intervention` is `NULL`, only the natural course is evaluated. A
  `natural` element is also added automatically when `ref_int` asks for
  the natural course and the list contains no `NULL` element; see
  `ref_int` and Details.

- ref_int:

  Reference intervention for contrasts. Either an integer index (`0` =
  natural course; `1`, `2`, … = elements of `intervention`) or a
  character name matching an element (e.g., `"always"`). `0` and
  `"natural"` both resolve to the `NULL` element of `intervention` if
  there is one, and otherwise add a `natural` element — which requires
  an exposure model in `models`. See Details. Default: `0`.

- time_len:

  length of the time in the data.
