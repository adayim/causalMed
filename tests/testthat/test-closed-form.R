# Closed-form (analytic truth) validation — pre-release review item V-3.
#
# Point-treatment DGP with binary V, A, M, Y. V enters BOTH the mediator and the
# outcome model, so the natural and the interventional cross-world functionals
# have *different* exact truths:
#
#   natural:        M ~ p(M | a', V)                 -- subject's own V
#   interventional: M ~ sum_v p(M | a', v) p(v)      -- V-marginal draw
#
# Both truths are computed by exact enumeration below (during the pre-release
# review they were additionally confirmed against a 5e6-row brute-force
# counterfactual simulation, agreeing to 3e-4). Each estimator is then checked
# against ITS OWN truth: this is the only test class in the suite where the
# target is exact and estimand-specific, so it catches an estimator silently
# computing the wrong functional.

expit <- stats::plogis
aV <- 0.5
pA <- function(V) expit(-0.2 + 0.4 * V)
pM <- function(A, V) expit(-0.5 + 1.2 * A + 0.3 * V)
pY <- function(A, M, V) expit(-1.0 + 0.8 * A + 1.0 * M + 0.5 * A * M + 0.3 * V)

grid_V <- c(0, 1)
pV     <- c(1 - aV, aV)

# E[Y_{a, M(a')}] with the mediator drawn conditionally on the subject's own V
phi_natural <- function(a, ap) {
  sum(vapply(seq_along(grid_V), function(i) {
    v <- grid_V[i]
    pV[i] * sum(vapply(c(0, 1), function(m)
      pY(a, m, v) * stats::dbinom(m, 1, pM(ap, v)), numeric(1)))
  }, numeric(1)))
}

# E[Y_{a, G(a')}] with the mediator drawn from its V-marginal distribution
phi_interventional <- function(a, ap) {
  pm1 <- sum(pV * pM(ap, grid_V))
  sum(vapply(seq_along(grid_V), function(i) {
    v <- grid_V[i]
    pV[i] * (pY(a, 1, v) * pm1 + pY(a, 0, v) * (1 - pm1))
  }, numeric(1)))
}

gen_cf_data <- function(n) {
  V <- stats::rbinom(n, 1, aV)
  A <- stats::rbinom(n, 1, pA(V))
  M <- stats::rbinom(n, 1, pM(A, V))
  Y <- stats::rbinom(n, 1, pY(A, M, V))
  data.frame(id = seq_len(n), time = 1, V = V, A = A, M = M, Y = Y)
}

cf_models <- function() {
  list(
    spec_model(A ~ V,               var_type = "binary", mod_type = "exposure"),
    spec_model(M ~ A + V,           var_type = "binary", mod_type = "mediator"),
    spec_model(Y ~ A + M + A:M + V, var_type = "binary", mod_type = "outcome")
  )
}

# Mean estimate over `nrep` independent datasets, plus its Monte Carlo SE.
replicate_estimates <- function(nrep, n, ..., .seed) {
  set.seed(.seed)
  datasets <- lapply(seq_len(nrep), function(i) gen_cf_data(n))
  est <- lapply(seq_along(datasets), function(i) {
    fit <- suppressWarnings(mediation(
      data = datasets[[i]], id_var = "id", base_vars = "V", exposure = "A",
      outcome = "Y", time_var = "time", models = cf_models(),
      R = 1, quiet = TRUE, seed = 1000 + i, ...
    ))
    stats::setNames(fit$effect_size$Est, fit$effect_size$Intervention)
  })
  m <- do.call(rbind, est)
  list(mean = colMeans(m),
       mcse = apply(m, 2, stats::sd) / sqrt(nrow(m)))
}

# Tolerance: the estimators are plug-ins, so they carry the ordinary
# finite-sample bias of a correctly-specified substitution estimator (measured
# at ~-0.003 for this DGP during the review). Allow 3 MC-SE plus that slack.
expect_close_to_truth <- function(res, truth, slack = 0.008) {
  for (nm in names(truth)) {
    tol <- 3 * res$mcse[[nm]] + slack
    testthat::expect_lt(abs(res$mean[[nm]] - truth[[nm]]), tol)
  }
}


testthat::test_that("the analytic truths are internally consistent", {
  # Natural and interventional cross-world functionals must NOT coincide here
  # (V confounds the mediator), otherwise the test could not tell them apart.
  testthat::expect_false(isTRUE(all.equal(phi_natural(1, 0),
                                          phi_interventional(1, 0))))
  # The a = a' functionals coincide by construction only in expectation over V;
  # each must be a valid probability.
  for (f in list(phi_natural, phi_interventional))
    for (a in c(0, 1)) for (ap in c(0, 1))
      testthat::expect_true(f(a, ap) > 0 && f(a, ap) < 1)
})


testthat::test_that("gcomp natural effects hit the exact natural truth", {
  testthat::skip_on_cran()
  truth <- c(Phi00 = phi_natural(0, 0), Phi10 = phi_natural(1, 0),
             Phi01 = phi_natural(0, 1), Phi11 = phi_natural(1, 1))
  res <- replicate_estimates(nrep = 15, n = 4000, .seed = 20260712,
                             mediation_type = "N", mc_sample = 20000L)
  expect_close_to_truth(res, truth)
})


testthat::test_that("gcomp interventional effects hit the exact interventional truth", {
  testthat::skip_on_cran()
  truth <- c(Phi00 = phi_interventional(0, 0),
             Phi10 = phi_interventional(1, 0),
             Phi11 = phi_interventional(1, 1))
  res <- replicate_estimates(nrep = 15, n = 4000, .seed = 20260713,
                             mediation_type = "I", mc_sample = 20000L, n_vw = 2L)
  expect_close_to_truth(res, truth)

  # The plug-in total effect uses the natural-course interventions, so it must
  # track the NATURAL contrast, not the interventional overall effect.
  te_truth <- phi_natural(1, 1) - phi_natural(0, 0)
  te_est   <- res$mean[["nat1"]] - res$mean[["nat0"]]
  testthat::expect_lt(abs(te_est - te_truth), 0.01)
})


testthat::test_that("TMLE natural effects hit the exact natural truth", {
  testthat::skip_on_cran()
  truth <- c(Phi00 = phi_natural(0, 0), Phi10 = phi_natural(1, 0),
             Phi01 = phi_natural(0, 1), Phi11 = phi_natural(1, 1))
  res <- replicate_estimates(nrep = 15, n = 4000, .seed = 20260714,
                             mediation_type = "N", estimator = "tmle")
  expect_close_to_truth(res, truth)
})
