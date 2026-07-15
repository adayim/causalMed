# Regression tests for the `seed` contract (pre-release review finding 3.1).
#
# The bug: get_args_for() drops NULL-valued arguments, so `seed = NULL` never
# reached .run_interventions(), whose own default (12345L) then took over. Two
# consequences, both silent:
#   (1) every `seed = NULL` run was seeded identically -> repeated runs returned
#       the same estimates instead of independent draws;
#   (2) set.seed(12345) fired inside .run_interventions() with no save/restore,
#       so the caller's GLOBAL RNG state was left deterministic -- data simulated
#       by the user AFTER a package call repeated across loop iterations.

seed_test_data <- function(n = 600, seed = 5L) {
  set.seed(seed)
  V <- rbinom(n, 1, 0.5)
  A <- rbinom(n, 1, plogis(0.2 * V))
  M <- rbinom(n, 1, plogis(-0.3 + 0.8 * A))
  Y <- rbinom(n, 1, plogis(-1 + 0.5 * A + 0.7 * M))
  data.frame(id = seq_len(n), time = 1, V = V, A = A, M = M, Y = Y)
}

seed_test_models <- function() {
  list(
    spec_model(A ~ V,         var_type = "binary", mod_type = "exposure"),
    spec_model(M ~ A + V,     var_type = "binary", mod_type = "mediator"),
    spec_model(Y ~ A + M + V, var_type = "binary", mod_type = "outcome")
  )
}

run_med <- function(d, seed) {
  suppressWarnings(mediation(
    data = d, id_var = "id", base_vars = "V", exposure = "A", outcome = "Y",
    time_var = "time", models = seed_test_models(), mediation_type = "I",
    mc_sample = 2000L, R = 1, quiet = TRUE, seed = seed
  ))$effect_size$Est
}

run_gform <- function(d, seed) {
  suppressWarnings(gformula(
    data = d, id_var = "id", base_vars = "V", exposure = "A",
    time_var = "time",
    models = list(
      spec_model(A ~ V,         var_type = "binary", mod_type = "exposure"),
      spec_model(Y ~ A + V,     var_type = "binary", mod_type = "outcome")
    ),
    intervention = list(natural = NULL, always = 1),
    mc_sample = 2000L, R = 1, quiet = TRUE, seed = seed
  ))$effect_size$Est
}


testthat::test_that("seed = NULL gives independent runs (not silently re-seeded)", {
  d <- seed_test_data()
  testthat::expect_false(isTRUE(all.equal(run_med(d, NULL), run_med(d, NULL))))
  testthat::expect_false(isTRUE(all.equal(run_gform(d, NULL), run_gform(d, NULL))))
})


testthat::test_that("seed = NULL does not leave the global RNG state deterministic", {
  d <- seed_test_data()

  # The core of the bug: from two DIFFERENT pre-call RNG states, the post-call
  # state used to be identical, so a user's own rnorm()/rbinom() after the call
  # repeated across loop iterations.
  set.seed(111); invisible(run_med(d, NULL)); post1 <- .Random.seed
  set.seed(222); invisible(run_med(d, NULL)); post2 <- .Random.seed
  testthat::expect_false(identical(post1, post2))

  # The user-visible symptom, stated directly.
  set.seed(111); invisible(run_med(d, NULL)); x1 <- runif(5)
  set.seed(222); invisible(run_med(d, NULL)); x2 <- runif(5)
  testthat::expect_false(isTRUE(all.equal(x1, x2)))

  set.seed(111); invisible(run_gform(d, NULL)); g1 <- .Random.seed
  set.seed(222); invisible(run_gform(d, NULL)); g2 <- .Random.seed
  testthat::expect_false(identical(g1, g2))
})


testthat::test_that("an integer seed reproduces exactly and restores the caller's RNG", {
  d <- seed_test_data()

  testthat::expect_equal(run_med(d, 42L), run_med(d, 42L))
  testthat::expect_equal(run_gform(d, 42L), run_gform(d, 42L))

  # A seeded call must not disturb the ambient stream (save/restore on exit).
  set.seed(777); before <- .Random.seed
  invisible(run_med(d, 42L))
  testthat::expect_identical(.Random.seed, before)

  set.seed(777); before_g <- .Random.seed
  invisible(run_gform(d, 42L))
  testthat::expect_identical(.Random.seed, before_g)

  # ... so draws taken after a seeded call continue the caller's own stream.
  set.seed(777); y1 <- runif(3)
  set.seed(777); invisible(run_med(d, 42L)); y2 <- runif(3)
  testthat::expect_equal(y1, y2)
})
