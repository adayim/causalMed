# Tests for user-supplied exposure regimes in mediation() and the machinery
# shared with gformula() (validator, fixed time grid).

# ---- check_static_regime() --------------------------------------------------

testthat::test_that("check_static_regime accepts length 1 or time_len and coerces to numeric", {
  f <- causalMed:::check_static_regime
  testthat::expect_identical(f(1, 5, "x"), 1)
  testthat::expect_identical(f(TRUE, 5, "x"), 1)
  testthat::expect_identical(f(c(0, 1, 1, 0, 0), 5, "x"), c(0, 1, 1, 0, 0))
  testthat::expect_identical(f(c(FALSE, TRUE, TRUE, FALSE, FALSE), 5, "x"),
                             c(0, 1, 1, 0, 0))
  # values = NULL: any numeric (gformula dose-style interventions)
  testthat::expect_identical(f(2.5, 5, "x"), 2.5)
  # values supplied and satisfied: returned unchanged
  testthat::expect_identical(f(c(0, 1, 1), 3, "x", values = c(0, 1)), c(0, 1, 1))
})

testthat::test_that("check_static_regime rejects bad type, length and values, naming the argument", {
  f <- causalMed:::check_static_regime
  testthat::expect_error(f("1", 5, "exposure_regime"),
                         "`exposure_regime` must be a numeric or logical vector", fixed = TRUE)
  testthat::expect_error(f(c(0, 1), 5, "exposure_regime"),
                         "`exposure_regime` must have length 1 or 5 (one value per distinct time point); got length 2",
                         fixed = TRUE)
  # length 0 is not a regime
  testthat::expect_error(f(numeric(0), 5, "x"), "must have length 1 or 5", fixed = TRUE)
  testthat::expect_error(f(c(0, 2, 1, 0, 0), 5, "reference_regime", values = c(0, 1)),
                         "`reference_regime` must take values in {0, 1}; got {2}", fixed = TRUE)
  testthat::expect_error(f(c(0, NA, 1, 0, 0), 5, "reference_regime", values = c(0, 1)),
                         "`reference_regime` must not contain NA", fixed = TRUE)
  # values = NULL (gformula): NA is still rejected, logical NA included.
  testthat::expect_error(f(c(1, NA, 1), 3, "intervention$x"),
                         "`intervention$x` must not contain NA", fixed = TRUE)
  testthat::expect_error(f(NA, 3, "intervention$x"),
                         "`intervention$x` must not contain NA", fixed = TRUE)
})

testthat::test_that("gformula's intervention check uses the shared validator", {
  data("nonsurvivaldata", package = "causalMed")
  m1 <- spec_model(A ~ V, var_type = "binary", mod_type = "exposure")
  m3 <- spec_model(Y_bin ~ A + V, var_type = "binary", mod_type = "outcome")
  models <- list(m1, m3)
  # nonsurvivaldata has 5 distinct time points
  testthat::expect_error(
    gformula(data = nonsurvivaldata, id_var = "id", base_vars = "V", exposure = "A",
             time_var = "time", models = models,
             intervention = list(natural = NULL, bad = c(1, 1, 1)),
             mc_sample = 100, R = 1, quiet = TRUE),
    "`intervention$bad` must have length 1 or 5", fixed = TRUE)
  testthat::expect_error(
    gformula(data = nonsurvivaldata, id_var = "id", base_vars = "V", exposure = "A",
             time_var = "time", models = models,
             intervention = list(natural = NULL, bad = "1"),
             mc_sample = 100, R = 1, quiet = TRUE),
    "`intervention$bad` must be a numeric or logical vector", fixed = TRUE)
  # An NA used to pass this check and abort deep in the simulation.
  testthat::expect_error(
    gformula(data = nonsurvivaldata, id_var = "id", base_vars = "V", exposure = "A",
             time_var = "time", models = models,
             intervention = list(natural = NULL, bad = c(1, NA, 1, 1, 1)),
             mc_sample = 100, R = 1, quiet = TRUE),
    "`intervention$bad` must not contain NA", fixed = TRUE)
  # Still accepted: scalar, full-length, logical, dyn_int, dose value
  testthat::expect_no_error(suppressWarnings(
    gformula(data = nonsurvivaldata, id_var = "id", base_vars = "V", exposure = "A",
             time_var = "time", models = models,
             intervention = list(natural = NULL, s = 1, v = c(0, 1, 1, 0, 0), l = TRUE,
                                 d = dyn_int(as.numeric(A > 0)), dose = 2),
             mc_sample = 100, R = 1, quiet = TRUE)))
})

# ---- Fixed time grid --------------------------------------------------------

# Times {1, 3, 4, 6}: the grid is the set of distinct observed values (four
# steps), not the range 1..6. Per-row binary outcome so an "outcome" model fits
# on every row; adequate for grid tests.
make_gap_data <- function(n = 300, tt = c(1, 3, 4, 6)) {
  set.seed(31)
  d <- data.table::CJ(id = seq_len(n), time = tt)
  data.table::setorderv(d, c("id", "time"))
  d[, V := rep(stats::rnorm(n), each = length(tt))]
  d[, A := stats::rbinom(.N, 1, 0.5)]
  d[, L := V + 0.5 * A + stats::rnorm(.N, 0, 0.5)]
  d[, Y := stats::rbinom(.N, 1, stats::plogis(-1 + 0.5 * A + 0.3 * L))]
  d[]
}
gap_models <- function() list(
  spec_model(A ~ V + time,     var_type = "binary", mod_type = "exposure"),
  spec_model(L ~ V + A + time, var_type = "normal", mod_type = "covariate"),
  spec_model(Y ~ V + A + L,    var_type = "binary", mod_type = "outcome"))

testthat::test_that("the grid is the distinct observed times: length 4 accepted, 6 rejected, last step at t = 6", {
  d <- make_gap_data()
  testthat::expect_error(
    gformula(data = d, id_var = "id", base_vars = "V", exposure = "A", time_var = "time",
             models = gap_models(), intervention = list(natural = NULL, v = rep(1, 6)),
             mc_sample = 100, R = 1, quiet = TRUE),
    "length 1 or 4", fixed = TRUE)
  fit <- suppressWarnings(gformula(
    data = d, id_var = "id", base_vars = "V", exposure = "A", time_var = "time",
    models = gap_models(), intervention = list(natural = NULL, v = c(0, 1, 1, 0)),
    mc_sample = 200, R = 1, quiet = TRUE, seed = 3, return_data = TRUE))
  testthat::expect_identical(unique(fit$sim_data$time), 6)
})

testthat::test_that(".run_interventions simulates the caller's grid, not the grid of the data it is handed", {
  d <- make_gap_data()
  d_short <- d[time != 6]                      # a resample that lost time 6
  run <- function(time_seq) {
    causalMed:::.run_interventions(
      data = d_short, id_var = "id", base_vars = "V", time_var = "time",
      exposure = "A", models = gap_models(),
      intervention = list(v = c(0, 1, 1, 0)),
      mc_sample = 100, return_data = TRUE, seed = 5, time_seq = time_seq)
  }
  # Fallback (no grid supplied): derived from the data handed in -> stops at 4.
  testthat::expect_identical(unique(run(NULL)$gform.data$v$time), 4)
  # Grid supplied by the caller: all four steps, last at 6.
  testthat::expect_identical(unique(run(c(1, 3, 4, 6))$gform.data$v$time), 6)
})

testthat::test_that("time_seq reaches the bootstrap replicates (plumbing; not a grid-shortening test)", {
  # Every resample of make_gap_data() contains all four time values (CJ gives
  # every subject every time), so this cannot detect a shortened grid -- test
  # (b) above covers that deterministically at the right seam. What this does
  # catch is time_seq being dropped or mangled on the way to a replicate.
  d <- make_gap_data()
  fit <- suppressWarnings(gformula(
    data = d, id_var = "id", base_vars = "V", exposure = "A", time_var = "time",
    models = gap_models(), intervention = list(natural = NULL, v = c(0, 1, 1, 0)),
    mc_sample = 100, R = 3, quiet = TRUE, seed = 3))
  bi <- fit$boot_estimates$interventions
  testthat::expect_equal(nrow(bi), 3L * 2L)
  testthat::expect_true(all(is.finite(bi$Est)))
})

# ---- intervention_spec machinery -------------------------------------------

testthat::test_that("intervention_spec accepts regime vectors and rejects bad ones", {
  is <- causalMed:::intervention_spec
  s <- is(c(0, 1, 1, 0, 0), list(M = c(0, 0, 0, 0, 0)))
  testthat::expect_s3_class(s, "causalMed_intervention")
  testthat::expect_identical(s$treatment, c(0, 1, 1, 0, 0))
  testthat::expect_identical(s$mediator_overrides$M, c(0, 0, 0, 0, 0))
  testthat::expect_identical(is(1)$treatment, 1)                     # scalar still fine
  testthat::expect_error(is("1"), "numeric vector", fixed = TRUE)
  testthat::expect_error(is(numeric(0)), "numeric vector", fixed = TRUE)
  testthat::expect_error(is(1, list(M = "0")), "numeric regime vector", fixed = TRUE)
  testthat::expect_error(is(1, list(0)), "named list", fixed = TRUE)
  testthat::expect_error(is(c(0, 1, 1), list(M = c(0, 1))),
                         "regime lengths must all be 1 or the same length", fixed = TRUE)
})

testthat::test_that("regime_key and slice_intervention_spec", {
  rk <- causalMed:::regime_key
  testthat::expect_identical(rk(c(0, 1, 1, 0, 0)), "0,1,1,0,0")
  testthat::expect_identical(rk(1), "1")
  testthat::expect_false(identical(rk(c(1, 1)), rk(1)))   # keys are exact, so recycle first
  # Strengthen: pin the key, not just its inequality.
  testthat::expect_identical(rk(c(1, 1)), "1,1")

  s <- causalMed:::intervention_spec(c(0, 1, 1, 0, 0),
                                     list(M1 = c(0, 0, 0, 0, 0), M2 = 1))
  k3 <- causalMed:::slice_intervention_spec(s, 3L)
  testthat::expect_s3_class(k3, "causalMed_intervention")
  testthat::expect_identical(k3$treatment, 1)
  testthat::expect_identical(k3$mediator_overrides, list(M1 = 0, M2 = 1))
  k1 <- causalMed:::slice_intervention_spec(s, 1L)
  testthat::expect_identical(k1$treatment, 0)

  # The nat0/nat1 path: empty overrides, sliced at every step of those
  # interventions.
  n <- causalMed:::slice_intervention_spec(
    causalMed:::intervention_spec(c(0, 1, 1, 0, 0), list()), 4L)
  testthat::expect_identical(n$treatment, 0)
  testthat::expect_identical(n$mediator_overrides, list())

  # Out-of-range k must ERROR, not silently give NA: the design relies on
  # `[[` here, and a "simplification" to `[` would inject NA into the
  # exposure column instead.
  testthat::expect_error(
    causalMed:::slice_intervention_spec(causalMed:::intervention_spec(c(0, 1, 1)), 5L))
})

testthat::test_that("build_mediation_interventions puts the regimes in the right slots", {
  b <- causalMed:::build_mediation_interventions
  a  <- c(0, 1, 1, 0, 0); a0 <- c(0, 0, 0, 0, 0)

  iv <- b(c("M1", "M2"), "I", regime = a, ref_regime = a0)
  testthat::expect_identical(names(iv), c("nat0", "nat1", "Phi00", "Phi10", "Phi1_1", "Phi11"))
  testthat::expect_identical(iv$nat0$treatment, a0)
  testthat::expect_identical(iv$nat1$treatment, a)
  testthat::expect_identical(iv$Phi00$treatment, a0)
  testthat::expect_identical(iv$Phi00$mediator_overrides, list(M1 = a0, M2 = a0))
  testthat::expect_identical(iv$Phi10$treatment, a)
  testthat::expect_identical(iv$Phi10$mediator_overrides, list(M1 = a0, M2 = a0))
  testthat::expect_identical(iv$Phi1_1$mediator_overrides, list(M1 = a, M2 = a0))
  testthat::expect_identical(iv$Phi11$mediator_overrides, list(M1 = a, M2 = a))

  nv <- b("M", "N", regime = a, ref_regime = a0)
  testthat::expect_identical(names(nv), c("Phi00", "Phi11", "Phi10", "Phi01"))
  testthat::expect_identical(nv$Phi10$treatment, a)
  testthat::expect_identical(nv$Phi10$mediator_overrides, list(M = a0))
  testthat::expect_identical(nv$Phi01$treatment, a0)
  testthat::expect_identical(nv$Phi01$mediator_overrides, list(M = a))

  # Mixed: vector regime with a SCALAR reference. This is the shape a
  # defaulted mediation() call would produce if it failed to recycle, and the
  # keys must then differ from the recycled form -- which is why mediation()
  # recycles both regimes before calling this.
  mx <- b("M", "I", regime = a, ref_regime = 0)
  testthat::expect_identical(mx$Phi00$mediator_overrides, list(M = 0))
  testthat::expect_false(
    identical(causalMed:::regime_key(mx$Phi00$treatment),
              causalMed:::regime_key(a0)))

  # Defaults are the old hardcoded scalars.
  d <- b("M", "I")
  testthat::expect_identical(d$Phi10$treatment, 1)
  testthat::expect_identical(d$Phi10$mediator_overrides, list(M = 0))
})

# ---- Engine: pools keyed by regime, regime applied per step ---------------

testthat::test_that("mediator pools are collected under the regime, step by step", {
  # Two collect_pool runs from the same seed whose regimes first differ at
  # index 3 must have IDENTICAL pool slices at indices 1-2 (same exposure,
  # same RNG consumption) and DIFFERENT slices at 3-4 (L ~ A and M ~ A + L).
  set.seed(4242)
  n <- 600; tt <- 0:3
  d <- data.table::CJ(id = seq_len(n), time = tt)
  data.table::setorderv(d, c("id", "time"))
  d[, V := rep(stats::rnorm(n), each = length(tt))]
  d[, A := stats::rbinom(.N, 1, 0.5)]
  d[, L := V + 0.5 * A + stats::rnorm(.N, 0, 0.5)]
  d[, M := 10 * time + 2 * A + L + stats::rnorm(.N, 0, 1)]
  d[, Y := stats::rbinom(.N, 1, stats::plogis(-2 + 0.5 * A + 0.02 * M + 0.3 * L))]
  d[, lag1_M := data.table::shift(M, fill = 0), by = id]
  mods <- list(
    spec_model(L ~ V + A + time,          var_type = "normal", mod_type = "covariate"),
    spec_model(M ~ A + L + lag1_M + time, var_type = "normal", mod_type = "mediator"),
    spec_model(Y ~ V + A + M + L,         var_type = "binary", mod_type = "outcome"))
  fit_mods <- causalMed:::fit_spec_models(mods, d)
  base_dat <- unique(d[, c("id", "V"), with = FALSE])

  run_pool <- function(regime) {
    set.seed(7)
    mc <- data.table::as.data.table(base_dat[sample(nrow(base_dat), 800, replace = TRUE), ])
    causalMed:::simulate_intervention(
      data = mc, models = fit_mods, exposure = "A", time_var = "time",
      time_seq = tt, intervention = causalMed:::intervention_spec(regime),
      init_recode = recodes(lag1_M = 0), in_recode = recodes(lag1_M = M),
      mediation_type = "I", collect_pool = TRUE)$pool$M
  }
  p_never <- run_pool(c(0, 0, 0, 0))
  p_late  <- run_pool(c(0, 0, 1, 1))
  testthat::expect_identical(p_never[[1]], p_late[[1]])
  testthat::expect_identical(p_never[[2]], p_late[[2]])
  testthat::expect_false(isTRUE(all.equal(p_never[[3]], p_late[[3]])))
  testthat::expect_false(isTRUE(all.equal(p_never[[4]], p_late[[4]])))
  # The shift is in the direction of the A effect (2 direct + 0.5 through L).
  testthat::expect_gt(mean(p_late[[3]]) - mean(p_never[[3]]), 1)
})

testthat::test_that("pools are keyed by the OVERRIDE regime, not the treatment regime", {
  # Discriminator: `cross` and `full` share a treatment regime and differ only
  # in which regime their mediator pool comes from. If pools were keyed by
  # treatment, both would draw from the same pool and the two Phi values would
  # coincide -- while the default scalar case, and the baseline fixture, would
  # still pass. This is the test that separates those two implementations.
  set.seed(4242)
  n <- 600; tt <- 0:3
  d <- data.table::CJ(id = seq_len(n), time = tt)
  data.table::setorderv(d, c("id", "time"))
  d[, V := rep(stats::rnorm(n), each = length(tt))]
  d[, A := stats::rbinom(.N, 1, 0.5)]
  d[, L := V + 0.5 * A + stats::rnorm(.N, 0, 0.5)]
  d[, M := 10 * time + 2 * A + L + stats::rnorm(.N, 0, 1)]
  d[, Y := stats::rbinom(.N, 1, stats::plogis(-2 + 0.5 * A + 0.02 * M + 0.3 * L))]
  d[, lag1_M := data.table::shift(M, fill = 0), by = id]
  mods <- list(
    spec_model(L ~ V + A + time,          var_type = "normal", mod_type = "covariate"),
    spec_model(M ~ A + L + lag1_M + time, var_type = "normal", mod_type = "mediator"),
    spec_model(Y ~ V + A + M + L,         var_type = "binary", mod_type = "outcome"))

  a  <- c(0, 1, 1, 0)
  a0 <- c(0, 0, 0, 0)
  res <- causalMed:::.run_interventions(
    data = d, id_var = "id", base_vars = "V", time_var = "time", exposure = "A",
    models = mods,
    intervention = list(
      ref   = causalMed:::intervention_spec(a0, list(M = a0)),
      cross = causalMed:::intervention_spec(a,  list(M = a0)),
      full  = causalMed:::intervention_spec(a,  list(M = a))),
    init_recode = recodes(lag1_M = 0), in_recode = recodes(lag1_M = M),
    mediation_type = "I", mc_sample = 800, n_vw = 1L, seed = 13, time_seq = tt)
  phi <- unlist(res$gform.data)
  testthat::expect_length(phi, 3L)
  testthat::expect_true(all(is.finite(phi)))
  # cross vs full differ ONLY in the mediator pool's regime.
  testthat::expect_false(isTRUE(all.equal(phi[["cross"]], phi[["full"]])))
  testthat::expect_false(isTRUE(all.equal(phi[["ref"]],   phi[["cross"]])))
})


# ---- mediation(): regimes end to end --------------------------------------

regime_models <- function() list(
  spec_model(L2    ~ A + V + time,      var_type = "binary", mod_type = "covariate"),
  spec_model(L1    ~ A + V + L2 + time, var_type = "normal", mod_type = "mediator"),
  spec_model(Y_bin ~ A + L1 + L2 + V,   var_type = "binary", mod_type = "outcome"))

run_med <- function(mtype, ...) {
  data("nonsurvivaldata", package = "causalMed")
  suppressWarnings(mediation(
    data = nonsurvivaldata, id_var = "id", base_vars = "V", exposure = "A",
    outcome = "Y_bin", time_var = "time", models = regime_models(),
    init_recode = recodes(lag1_A = 0), in_recode = recodes(lag1_A = A),
    mediation_type = mtype, mc_sample = 400L, quiet = TRUE, seed = 11L, ...))
}

testthat::test_that("a mixed regime runs under I and N and is recorded recycled in all.args", {
  for (mtype in c("I", "N")) {
    fit <- run_med(mtype, R = 1L, exposure_regime = c(0, 1, 1, 0, 0), reference_regime = 0)
    testthat::expect_identical(fit$all.args$exposure_regime,  c(0, 1, 1, 0, 0))
    testthat::expect_identical(fit$all.args$reference_regime, c(0, 0, 0, 0, 0))
    testthat::expect_true(all(is.finite(fit$effect_size$Est)))
    testthat::expect_true(all(c("Direct effect", "Indirect effect", "Total effect") %in%
                                fit$estimate$Effect))
  }
})

testthat::test_that("the regime is applied step by step (final-step snapshot)", {
  # Regimes whose LAST element differs from their first, so the end-of-follow-up
  # snapshot separates per-step application from any scalar fallback.
  a  <- c(1, 1, 1, 1, 0)
  a0 <- c(0, 0, 0, 0, 1)
  for (mtype in c("I", "N")) {
    fit <- run_med(mtype, R = 1L, n_vw = 1L, return_data = TRUE,
                   exposure_regime = a, reference_regime = a0)
    snap <- fit$sim_data[, .(A = unique(A), lag1_A = unique(lag1_A)), by = Intervention]
    testthat::expect_true(all(vapply(split(snap, snap$Intervention), nrow, integer(1)) == 1L),
                          label = "one distinct (A, lag1_A) pair per intervention")
    on_a  <- if (mtype == "I") c("nat1", "Phi10", "Phi11") else c("Phi11", "Phi10")
    on_a0 <- if (mtype == "I") c("nat0", "Phi00")          else c("Phi00", "Phi01")
    # Pin the names first: `snap[Intervention %in% <typo>, A]` is numeric(0) and
    # `all(numeric(0) == x)` is TRUE, so a renamed intervention would silently
    # drop the assertions below rather than fail.
    testthat::expect_setequal(snap$Intervention, c(on_a, on_a0))
    testthat::expect_true(all(snap[Intervention %in% on_a,  A] == a[5]))
    testthat::expect_true(all(snap[Intervention %in% on_a0, A] == a0[5]))
    # lag1_A at the final step is A at step 4
    testthat::expect_true(all(snap[Intervention %in% on_a,  lag1_A] == a[4]))
    testthat::expect_true(all(snap[Intervention %in% on_a0, lag1_A] == a0[4]))
  }
})

testthat::test_that("regimes reach the bootstrap replicates and two mediators work", {
  fit <- run_med("I", R = 2L, exposure_regime = c(0, 1, 1, 0, 0))
  testthat::expect_equal(nrow(fit$boot_estimates$interventions), 2L * 5L)
  testthat::expect_true(all(is.finite(fit$boot_estimates$interventions$Est)))

  data("nonsurvivaldata", package = "causalMed")
  mods2 <- list(
    spec_model(L2    ~ A + V + time,           var_type = "binary", mod_type = "covariate"),
    spec_model(L1    ~ A + V + L2 + time,      var_type = "normal", mod_type = "mediator"),
    spec_model(M     ~ A + V + L1 + L2 + time, var_type = "normal", mod_type = "mediator"),
    spec_model(Y_bin ~ A + L1 + M + L2 + V,    var_type = "binary", mod_type = "outcome"))
  fit2 <- suppressWarnings(mediation(
    data = nonsurvivaldata, id_var = "id", base_vars = "V", exposure = "A",
    outcome = "Y_bin", time_var = "time", models = mods2,
    mediation_type = "I", exposure_regime = c(0, 1, 1, 0, 0),
    mc_sample = 300L, R = 1L, quiet = TRUE, seed = 11L))
  testthat::expect_identical(fit2$effect_size$Intervention,
                             c("nat0", "nat1", "Phi00", "Phi10", "Phi1_1", "Phi11"))
  testthat::expect_true(all(is.finite(fit2$effect_size$Est)))
})

testthat::test_that("a mediator with no collected pool errors instead of silently changing estimand", {
  # A pool-key miss used to fall through to simulate_data()'s within-step
  # permutation -- a different estimand, with no error.
  set.seed(99)
  tt <- 0:2
  d <- data.table::CJ(id = seq_len(200), time = tt)
  data.table::setorderv(d, c("id", "time"))
  d[, V := rep(stats::rnorm(200), each = length(tt))]
  d[, A := stats::rbinom(.N, 1, 0.5)]
  d[, M := A + V + stats::rnorm(.N)]
  d[, Y := stats::rbinom(.N, 1, stats::plogis(-1 + 0.3 * M))]
  mods <- list(
    spec_model(M ~ A + V,     var_type = "normal", mod_type = "mediator"),
    spec_model(Y ~ A + M + V, var_type = "binary", mod_type = "outcome"))
  testthat::expect_error(
    causalMed:::.run_interventions(
      data = d, id_var = "id", base_vars = "V", time_var = "time", exposure = "A",
      models = mods,
      # Names a variable that is not a fitted mediator response: a pool is
      # collected under the regime, but it holds no slice for "Z".
      intervention = list(bad = causalMed:::intervention_spec(c(1, 1, 1),
                                                             list(Z = c(0, 0, 0)))),
      mediation_type = "I", mc_sample = 200, n_vw = 1L, seed = 1, time_seq = tt),
    "no mediator pool for", fixed = TRUE)

  # One layer down: simulate_data() used to fall back to a within-step
  # permutation when handed an override with no pool slice. It now stops.
  fit_mods <- causalMed:::fit_spec_models(mods, d)
  mc <- data.table::as.data.table(unique(d[, c("id", "V"), with = FALSE]))
  testthat::expect_error(
    causalMed:::simulate_intervention(
      data = mc, models = fit_mods, exposure = "A", time_var = "time",
      time_seq = tt,
      intervention = causalMed:::intervention_spec(c(1, 1, 1), list(M = c(0, 0, 0))),
      mediation_type = "I", med_pool = NULL),
    "no mediator pool slice for 'M'", fixed = TRUE)
})

# ---- Support diagnostic -----------------------------------------------------

testthat::test_that("regime_support counts subjects whose observed exposure equals the regime at every observed time", {
  d <- data.frame(
    id   = c(1, 1, 1,  2, 2, 2,  3, 3,  4, 4, 4,  5, 5, 5),
    time = c(0, 1, 2,  0, 1, 2,  0, 1,  0, 1, 2,  0, 2, 1),   # id 5 rows out of order
    A    = c(0, 1, 1,  0, 0, 0,  0, 1,  0, NA, 1,  1, 1, 0))
  rs <- causalMed:::regime_support(d, "id", "time", "A", time_seq = 0:2,
                                   regimes = list("a" = c(0, 1, 1), "a*" = c(0, 0, 0)))
  testthat::expect_identical(rs$regime, c("a", "a*"))
  testthat::expect_identical(rs$values, c("0 1 1", "0 0 0"))
  # a : id 1 exactly; id 3 on its two observed times (short follow-up); id 4 has NA; id 5 no.
  # a*: id 2 only.
  testthat::expect_identical(rs$n_following, c(2L, 1L))
  testthat::expect_equal(rs$prop_following, c(2 / 5, 1 / 5))
  # n_complete: of those following, only subjects also observed at all three
  # grid times (0, 1, 2) count. id 1 (a) and id 2 (a*) are complete; id 3 (a)
  # is only observed at two of the three grid times, so it does not count.
  testthat::expect_identical(rs$n_complete, c(1L, 1L))
})

testthat::test_that("regime_support counts distinct grid times and ignores NA-time rows", {
  d <- data.frame(
    id   = c(1, 1, 1, 1,   2, 2, 2, 2, 2,   3, 3, 3, 3),
    time = c(0, 0, 1, 2,   0, 1, 2, 3, NA,  0, 1, 2, 3),
    A    = c(1, 1, 1, 1,   1, 1, 1, 1, 0,   0, 0, 0, 0))
  rs <- causalMed:::regime_support(d, "id", "time", "A", time_seq = 0:3,
                                   regimes = list("a" = rep(1, 4), "a*" = rep(0, 4)))
  # id 1: four rows but only three distinct times (no t = 3) -> follows a,
  #       NOT complete (a row count would call it complete).
  # id 2: its NA-time row (A = 0) has no grid position and is ignored, so it
  #       follows a and is complete.
  # id 3: follows a*, complete.
  testthat::expect_identical(rs$n_following, c(2L, 1L))
  testthat::expect_identical(rs$n_complete,  c(1L, 1L))
  testthat::expect_equal(rs$prop_following, c(2 / 3, 1 / 3))
})

testthat::test_that("mediation attaches regime_support and time_seq to data_summary", {
  fit <- run_med("I", R = 1L, exposure_regime = c(0, 1, 1, 0, 0))
  ds <- fit$data_summary
  testthat::expect_identical(ds$time_seq, c(0, 1, 2, 3, 4))
  rs <- ds$regime_support
  testthat::expect_s3_class(rs, "data.frame")
  testthat::expect_identical(names(rs), c("regime", "values", "n_following",
                                          "prop_following", "n_complete"))
  testthat::expect_identical(rs$values, c("0 1 1 0 0", "0 0 0 0 0"))
  testthat::expect_true(all(rs$n_following >= 0L))
  testthat::expect_equal(rs$prop_following, rs$n_following / ds$n_id)
})

testthat::test_that("n_complete equals n_following on a balanced panel and is smaller under censoring", {
  # nonsurvivaldata: every subject has all 5 times, so following at every
  # observed time IS following at every grid time -- the two counts coincide.
  data("nonsurvivaldata", package = "causalMed")
  rs_bal <- causalMed:::regime_support(
    nonsurvivaldata, "id", "time", "A", time_seq = sort(unique(nonsurvivaldata$time)),
    regimes = list("a" = c(0, 1, 1, 0, 0)))
  testthat::expect_identical(rs_bal$n_complete, rs_bal$n_following)

  # survivaldata: heavy right-censoring, so most subjects who match on their
  # (short) observed history are not observed for the whole grid.
  data("survivaldata", package = "causalMed")
  rs_cens <- causalMed:::regime_support(
    survivaldata, "id", "time", "A", time_seq = sort(unique(survivaldata$time)),
    regimes = list("a" = c(0, 1, 1, 0, 0)))
  testthat::expect_true(rs_cens$n_complete < rs_cens$n_following)
})


# ---- print() ----------------------------------------------------------------

testthat::test_that("print shows the regimes and support lines; defaults keep the old legend", {
  fit_mix <- run_med("I", R = 1L, exposure_regime = c(0, 1, 1, 0, 0))
  out <- paste(utils::capture.output(print(fit_mix)), collapse = "\n")
  testthat::expect_match(out, "Exposure regime  a : 0 1 1 0 0   (time = 0 1 2 3 4)", fixed = TRUE)
  testthat::expect_match(out, "Reference regime a*: 0 0 0 0 0", fixed = TRUE)
  testthat::expect_match(out, "Observed subjects following a  (0 1 1 0 0):", fixed = TRUE)
  testthat::expect_match(out, "Observed subjects following a* (0 0 0 0 0):", fixed = TRUE)
  testthat::expect_match(out, "Phi10 = E[Y(a, G_a*)]:  exposure = a, mediators ~ a* pool  [cross-regime]", fixed = TRUE)
  testthat::expect_false(grepl("exposure=1", out, fixed = TRUE))

  fit_def <- run_med("I", R = 1L)
  out_def <- paste(utils::capture.output(print(fit_def)), collapse = "\n")
  testthat::expect_false(grepl("Exposure regime", out_def, fixed = TRUE))
  testthat::expect_match(out_def, "Observed subjects following a  (1 1 1 1 1):", fixed = TRUE)
  # The whole default legend, byte for byte: this text is what README.md and
  # the vignettes capture, so any reword must be a deliberate, visible change.
  testthat::expect_match(out_def, paste0(
    "  Under interventional effects, each intervention draws its mediators from independently-permuted pools (G):\n",
    "  Phi11 = E[Y(a=1, G1)]:  exposure=1, mediators ~ a=1 pool  [reference]\n",
    "  Phi10 = E[Y(a=1, G0)]:  exposure=1, mediators ~ a=0 pool  [cross-regime]\n",
    "  Phi00 = E[Y(a=0, G0)]:  exposure=0, mediators ~ a=0 pool  [reference]\n",
    "  nat1/nat0 = E[Y(a=1)]/E[Y(a=0)]:  exposure fixed, mediators natural (used for the total effect)"),
    fixed = TRUE)

  fit_n <- run_med("N", R = 1L, exposure_regime = c(0, 1, 1, 0, 0))
  out_n <- paste(utils::capture.output(print(fit_n)), collapse = "\n")
  testthat::expect_match(out_n, "Phi10 = E[Y(a, M(a*))]:  exposure = a, mediator under a*  [cross-world]", fixed = TRUE)
  # The decomposition legend follows the regimes too (it used to stay
  # hardcoded to always vs never under "N").
  testthat::expect_match(out_n,
    "  Total effect    = Phi11 - Phi00 =  E[Y(a,M(a))] - E[Y(a*,M(a*))]\n",
    fixed = TRUE)
  testthat::expect_match(out_n,
    "  Direct effect   = Phi10 - Phi00 =  E[Y(a,M(a*))] - E[Y(a*,M(a*))]\n",
    fixed = TRUE)
  testthat::expect_match(out_n,
    "  Indirect effect = Phi11 - Phi10 =  E[Y(a,M(a))] - E[Y(a,M(a*))]\n",
    fixed = TRUE)
  testthat::expect_false(grepl("E[Y(1,M(1))]", out_n, fixed = TRUE))

  # Default "N": both legends byte for byte as before.
  out_nd <- paste(utils::capture.output(print(run_med("N", R = 1L))), collapse = "\n")
  testthat::expect_match(out_nd, paste0(
    "  Phi11 = E[Y(a=1, M(1))]:  exposure=1, mediator under a=1\n",
    "  Phi10 = E[Y(a=1, M(0))]:  exposure=1, mediator under a=0  [cross-world]\n",
    "  Phi01 = E[Y(a=0, M(1))]:  exposure=0, mediator under a=1  [cross-world]\n",
    "  Phi00 = E[Y(a=0, M(0))]:  exposure=0, mediator under a=0\n"), fixed = TRUE)
  # Phi01 is simulated and printed, so the legend must describe it.
  testthat::expect_match(out_nd, "Phi01 is reported for completeness", fixed = TRUE)
  testthat::expect_match(out_nd, paste0(
    "  Total effect    = Phi11 - Phi00 =  E[Y(1,M(1))] - E[Y(0,M(0))]\n",
    "  Direct effect   = Phi10 - Phi00 =  E[Y(1,M(0))] - E[Y(0,M(0))]\n",
    "  Indirect effect = Phi11 - Phi10 =  E[Y(1,M(1))] - E[Y(1,M(0))]\n"), fixed = TRUE)

  # Two mediators: the sequential Phi1_k line appears in both legend branches.
  data("nonsurvivaldata", package = "causalMed")
  mods2 <- list(
    spec_model(L2    ~ A + V + time,           var_type = "binary", mod_type = "covariate"),
    spec_model(L1    ~ A + V + L2 + time,      var_type = "normal", mod_type = "mediator"),
    spec_model(M     ~ A + V + L1 + L2 + time, var_type = "normal", mod_type = "mediator"),
    spec_model(Y_bin ~ A + L1 + M + L2 + V,    var_type = "binary", mod_type = "outcome"))
  med2 <- function(...) suppressWarnings(mediation(
    data = nonsurvivaldata, id_var = "id", base_vars = "V", exposure = "A",
    outcome = "Y_bin", time_var = "time", models = mods2, mediation_type = "I",
    mc_sample = 300L, R = 1L, quiet = TRUE, seed = 11L, ...))
  o2d <- paste(utils::capture.output(print(med2())), collapse = "\n")
  testthat::expect_match(o2d,
    "  Phi1_k:  exposure=1, first k mediators ~ a=1 pool, rest ~ a=0 pool  [sequential]",
    fixed = TRUE)
  o2m <- paste(utils::capture.output(
    print(med2(exposure_regime = c(0, 1, 1, 0, 0)))), collapse = "\n")
  testthat::expect_match(o2m,
    "  Phi1_k:  exposure = a, first k mediators ~ a pool, rest ~ a* pool  [sequential]",
    fixed = TRUE)
})

testthat::test_that("print survives objects saved before regimes existed", {
  fit <- run_med("I", R = 1L)
  fit$all.args$exposure_regime  <- NULL
  fit$all.args$reference_regime <- NULL
  fit$data_summary$regime_support <- NULL
  fit$data_summary$time_seq <- NULL
  testthat::expect_no_error(utils::capture.output(print(fit)))
})

testthat::test_that("long regimes and time grids are elided in print, full in the object", {
  s <- paste(c(rep(0, 40), rep(1, 60)), collapse = " ")
  testthat::expect_identical(causalMed:::.fmt_regime_str(s), "0 x 40, 1 x 60")
  # Short regimes are untouched.
  testthat::expect_identical(causalMed:::.fmt_regime_str("0 1 1 0 0"), "0 1 1 0 0")
  # A regime that switches often: RLE would be LONGER than the raw string, so
  # it is elided head/tail instead -- and never printed longer than it was.
  alt <- paste(rep(0:1, length.out = 21L), collapse = " ")
  alt_fmt <- causalMed:::.fmt_regime_str(alt)
  testthat::expect_lt(nchar(alt_fmt), nchar(alt))
  testthat::expect_identical(alt_fmt, "0 1 0 ... 0 1 0  (21 points)")
  testthat::expect_identical(causalMed:::.fmt_times(0:4), "0 1 2 3 4")
  testthat::expect_match(causalMed:::.fmt_times(0:1824), "(1825 points)", fixed = TRUE)
})

# ---- Guards that must fail closed -------------------------------------------

testthat::test_that("default_regime_pair() needs both regimes, not an empty one", {
  f <- causalMed:::default_regime_pair
  # Legacy objects, saved before the arguments existed, carry neither.
  testthat::expect_true(f(NULL, NULL))
  testthat::expect_true(f(c(1, 1, 1), c(0, 0, 0)))
  testthat::expect_false(f(c(1, 0, 1), c(0, 0, 0)))
  # all(logical(0)) is TRUE, so a half-populated pair must not read as the
  # default: it gates estimator = "tmle" and picks the print legend.
  testthat::expect_false(f(c(1, 1, 1), NULL))
  testthat::expect_false(f(c(1, 1, 1), numeric(0)))
  testthat::expect_false(f(numeric(0), c(0, 0, 0)))
})

testthat::test_that("check_static_regime() rejects a matrix instead of flattening it", {
  f <- causalMed:::check_static_regime
  testthat::expect_error(
    f(matrix(c(0, 1, 1), nrow = 1), 3, "exposure_regime", values = c(0, 1)),
    "`exposure_regime` must be a plain vector, not a matrix", fixed = TRUE)
  testthat::expect_error(
    f(array(c(0, 1, 1, 0), dim = c(2, 1, 2)), 4, "reference_regime"),
    "must be a plain vector, not an array", fixed = TRUE)
  # A plain vector of the same values still passes.
  testthat::expect_identical(f(c(0, 1, 1), 3, "exposure_regime", values = c(0, 1)),
                             c(0, 1, 1))
})

testthat::test_that("regime_support() counts the same subjects in both parts of the proportion", {
  f  <- causalMed:::regime_support
  dt <- data.table::data.table(id = c("1", "1", NA, NA), time = c(0, 1, 0, 1),
                               A = c(1, 1, 1, 1))
  # tapply() discards the NA-id group, so counting it in the denominator would
  # report 1 of 2 for a regime both rows follow.
  out <- f(dt, "id", "time", "A", 0:1, list(a = c(1, 1)))
  testthat::expect_identical(out$n_following, 1L)
  testthat::expect_identical(out$prop_following, 1)
})

testthat::test_that("simulate_data() refuses an unsliced regime", {
  testthat::expect_error(
    causalMed:::simulate_data(
      data         = data.table::data.table(id = 1:5, A = 0, L = 0),
      exposure     = "A",
      models       = list(),
      intervention = causalMed:::intervention_spec(c(0, 1, 1), list())),
    "the exposure regime must be a single value, got 3", fixed = TRUE)
})

testthat::test_that("bootstrap_helper() rejects a missing time_seq at the call site", {
  testthat::expect_error(
    causalMed:::bootstrap_helper(
      data = data.table::data.table(id = 1L, time = 0L, A = 0L, Y = 0L),
      id_var = "id", base_vars = character(0), time_var = "time",
      exposure = "A", models = list(), intervention = NULL,
      R = 1L, progress_bar = FALSE),
    "`time_seq` is required", fixed = TRUE)
})
