# Tests for mediation(estimator = "tmle") — Zheng & van der Laan (2017) TMLE
# for longitudinal natural direct/indirect effects (mediation_type = "N").

# Small DGP with good positivity (coin-flip-ish exposure) so every regime has
# support: T = 3, binary L (post-mediator), binary M, binary end-of-follow-up Y.
gen_tmle_testdata <- function(n, seed = 1L) {
  set.seed(seed)
  expit <- plogis
  V <- rbinom(n, 1, 0.5)
  A_prev <- rep(0, n); M_prev <- rep(0, n); L_prev <- rep(0, n)
  out <- vector("list", 3L)
  for (t in 1:3) {
    A <- rbinom(n, 1, expit(0.2 * V + 0.3 * L_prev - 0.1 * A_prev))
    M <- rbinom(n, 1, expit(-0.4 + 0.8 * A + 0.3 * V + 0.3 * M_prev))
    L <- rbinom(n, 1, expit(-0.2 + 0.5 * A + 0.4 * M - 0.2 * L_prev))
    Y <- if (t == 3) rbinom(n, 1, expit(-1.2 + 0.6 * A + 0.7 * M + 0.4 * L +
                                          0.3 * V)) else rep(NA_real_, n)
    out[[t]] <- data.frame(id = seq_len(n), time = t, V = V,
                           A = A, M = M, L = L, Y = Y,
                           lag_A = A_prev, lag_M = M_prev, lag_L = L_prev)
    A_prev <- A; M_prev <- M; L_prev <- L
  }
  do.call(rbind, out)
}

tmle_test_models <- function() {
  list(
    spec_model(A ~ V + lag_A + lag_L + time, var_type = "binary",
               mod_type = "exposure"),
    spec_model(M ~ V + A + lag_M + time,     var_type = "binary",
               mod_type = "mediator"),
    spec_model(L ~ V + A + M + lag_L + time, var_type = "binary",
               mod_type = "covariate"),
    spec_model(Y ~ V + A + M + L,            var_type = "binary",
               mod_type = "outcome")
  )
}

tmle_test_args <- function(d) {
  list(data = d, id_var = "id", base_vars = "V",
       exposure = "A", outcome = "Y", time_var = "time",
       models = tmle_test_models(),
       init_recode = recodes(lag_A = 0, lag_M = 0, lag_L = 0),
       in_recode   = recodes(lag_A = A, lag_M = M, lag_L = L),
       mediation_type = "N", R = 1, quiet = TRUE)
}

testthat::test_that("tmle estimator returns a well-formed gformula object", {
  testthat::skip_on_cran()
  d <- gen_tmle_testdata(2000)
  fit <- suppressWarnings(
    do.call(mediation, c(tmle_test_args(d), list(estimator = "tmle")))
  )

  testthat::expect_s3_class(fit, "gformula")
  testthat::expect_setequal(fit$effect_size$Intervention,
                            c("Phi00", "Phi10", "Phi01", "Phi11"))
  testthat::expect_true(all(c("Sd", "norm_lcl", "norm_ucl") %in%
                              names(fit$effect_size)))
  testthat::expect_true(all(is.finite(fit$effect_size$Est)))
  testthat::expect_true(all(fit$effect_size$Est > 0 & fit$effect_size$Est < 1))
  testthat::expect_true(all(c("Indirect effect", "Direct effect",
                              "Total effect") %in% fit$estimate$Effect))
  # Natural-effect decomposition is exact: NDE + NIE = TE
  rd <- setNames(fit$estimate$RD, fit$estimate$Effect)
  testthat::expect_equal(rd[["Indirect effect"]] + rd[["Direct effect"]],
                         rd[["Total effect"]], tolerance = 1e-10)
  # EIC diagnostics present
  testthat::expect_named(fit$tmle_diag$eic_mean,
                         c("Phi00", "Phi10", "Phi01", "Phi11"))
  testthat::expect_output(print(fit), "TMLE")

  # The test DGP has L depend on A (intermediate confounder), so the
  # identifiability caveat must be retained on the object and re-surfaced by
  # print() in a short form.
  testthat::expect_true("L" %in% fit$intermediate_confounders)
  testthat::expect_output(print(fit), "Identifiability")
  testthat::expect_output(print(fit), "not point-identified|NOT point-identified")
})

testthat::test_that("tmle solves the EIC equation on well-supported data", {
  testthat::skip_on_cran()
  d <- gen_tmle_testdata(2000, seed = 7L)
  fit <- suppressWarnings(
    do.call(mediation, c(tmle_test_args(d), list(estimator = "tmle")))
  )
  # Post-targeting, the empirical mean of each EIC should be ~0 relative to
  # its sampling scale (sd(EIC)/sqrt(n)).
  eic_scale <- apply(fit$tmle_diag$eic, 2, sd) / sqrt(fit$tmle_diag$n)
  testthat::expect_true(all(abs(fit$tmle_diag$eic_mean) < 0.25 * eic_scale))
})

testthat::test_that("tmle and gcomp agree on well-specified, well-supported data", {
  testthat::skip_on_cran()
  d <- gen_tmle_testdata(4000, seed = 11L)
  fit_t <- suppressWarnings(
    do.call(mediation, c(tmle_test_args(d), list(estimator = "tmle")))
  )
  fit_g <- suppressWarnings(
    do.call(mediation, c(tmle_test_args(d),
                         list(estimator = "gcomp", mc_sample = 50000,
                              seed = 99L)))
  )
  pt <- setNames(fit_t$effect_size$Est, fit_t$effect_size$Intervention)
  pg <- setNames(fit_g$effect_size$Est, fit_g$effect_size$Intervention)
  sds <- setNames(fit_t$effect_size$Sd, fit_t$effect_size$Intervention)
  for (ph in names(pt)) {
    testthat::expect_lt(abs(pt[[ph]] - pg[[ph]]), 3 * sds[[ph]] + 0.01)
  }
})

testthat::test_that("tmle input validation errors are informative", {
  testthat::skip_on_cran()
  d <- gen_tmle_testdata(200)
  args <- tmle_test_args(d)

  # interventional path has no TMLE
  a1 <- args; a1$mediation_type <- "I"
  testthat::expect_error(
    do.call(mediation, c(a1, list(estimator = "tmle"))),
    "only available for mediation_type = 'N'"
  )

  # exposure model required
  a2 <- args; a2$models <- args$models[-1]
  testthat::expect_error(
    suppressWarnings(do.call(mediation, c(a2, list(estimator = "tmle")))),
    "requires an exposure model"
  )

  # subset not supported
  a3 <- args
  a3$models[[3]] <- spec_model(L ~ V + A + M + lag_L + time,
                               var_type = "binary", mod_type = "covariate",
                               subset = time > 1)
  testthat::expect_error(
    suppressWarnings(do.call(mediation, c(a3, list(estimator = "tmle")))),
    "does not support spec_model\\(subset"
  )

  # non-binary outcome rejected
  a4 <- args
  a4$data$Y <- ifelse(is.na(a4$data$Y), NA, a4$data$Y + 0.5)
  testthat::expect_error(
    suppressWarnings(do.call(mediation, c(a4, list(estimator = "tmle")))),
    "binary \\{0, 1\\} outcome"
  )
})
