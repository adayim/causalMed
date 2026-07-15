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


testthat::test_that("tmle rejects recodes it cannot honour (review finding 3.2)", {
  testthat::skip_on_cran()
  d <- gen_tmle_testdata(200)
  d$day <- d$time
  args <- tmle_test_args(d)

  # out_recode is never evaluated by the targeted engine -- must not be
  # silently ignored.
  a1 <- args
  a1$out_recode <- recodes(lag_A = A)
  testthat::expect_error(
    suppressWarnings(do.call(mediation, c(a1, list(estimator = "tmle")))),
    "does not support out_recode"
  )

  # Derived (non-lag) in_recode entries would be dropped from the lag map,
  # leaving stale values in the R-step and un-intervened exposure-derived
  # columns in the clever covariates.
  a2 <- args
  a2$in_recode <- recodes(lag_A = A, lag_M = M, lag_L = L, daysq = day^2)
  testthat::expect_error(
    suppressWarnings(do.call(mediation, c(a2, list(estimator = "tmle")))),
    "only supports in_recode entries that copy a single column"
  )

  # init_recode must be a constant or a bare column name.
  a3 <- args
  a3$init_recode <- recodes(lag_A = 0, lag_M = 0, lag_L = 0, daysq = day^2)
  testthat::expect_error(
    suppressWarnings(do.call(mediation, c(a3, list(estimator = "tmle")))),
    "only supports init_recode entries"
  )

  # Chained exposure lags (lag2_A = lag_A) are bare names, so they pass the
  # single-column check, but .tmle_set_regime intervenes only first-order
  # exposure lags -- they must be rejected, not silently left observed
  # (review 2026-07-15 finding 3a).
  a4 <- args
  a4$data$lag2_A <- 0
  a4$in_recode   <- recodes(lag_A = A, lag_M = M, lag_L = L, lag2_A = lag_A)
  a4$init_recode <- recodes(lag_A = 0, lag_M = 0, lag_L = 0, lag2_A = 0)
  testthat::expect_error(
    suppressWarnings(do.call(mediation, c(a4, list(estimator = "tmle")))),
    "first-order exposure lags"
  )
  # ... while chains of NON-exposure columns remain allowed (covariate
  # history is conditioned on as observed; the lag-map shift handles depth).
  a5 <- args
  a5$data$lag2_M <- 0
  a5$in_recode   <- recodes(lag_A = A, lag_M = M, lag_L = L, lag2_M = lag_M)
  a5$init_recode <- recodes(lag_A = 0, lag_M = 0, lag_L = 0, lag2_M = 0)
  testthat::expect_no_error(
    suppressWarnings(do.call(mediation, c(a5, list(estimator = "tmle"))))
  )

  # The lag-style recodes actually used by the TMLE examples still pass.
  testthat::expect_no_error(
    suppressWarnings(do.call(mediation, c(args, list(estimator = "tmle"))))
  )
})


testthat::test_that("tmle propagation intervenes exposure lags at t=2 (review 2026-07-15 finding 1)", {
  testthat::skip_on_cran()

  # T = 2 DGP whose R-block covariate loads directly on lag_A (coef 1.4): the
  # R-step marginalization then genuinely depends on exposure history -- the
  # one structural family the lag-Markov DGPs above cannot see. The buggy
  # propagation (exposure lag left at its init value instead of the regime
  # value) biased Phi10/Phi11 by ~0.04 (~7 SE at this n) on this DGP. Truth
  # is exact: only the t = 2 nodes feed Y.
  expit <- stats::plogis
  gen_lag_data <- function(n, seed) {
    set.seed(seed)
    V  <- rbinom(n, 1, 0.5)
    A1 <- rbinom(n, 1, expit(-0.2 + 0.4 * V))
    R1 <- rbinom(n, 1, expit(-0.4 + 0.9 * A1 + 0.3 * V))
    Z1 <- rbinom(n, 1, expit(-0.3 + 0.8 * A1 + 0.7 * R1 + 0.3 * V))
    A2 <- rbinom(n, 1, expit(-0.2 + 0.4 * V + 0.6 * A1))
    R2 <- rbinom(n, 1, expit(-0.4 + 0.9 * A2 + 1.4 * A1 + 0.3 * V))
    Z2 <- rbinom(n, 1, expit(-0.3 + 0.8 * A2 + 0.7 * R2 + 0.3 * V))
    Y2 <- rbinom(n, 1, expit(-0.5 + 0.7 * A2 + 0.9 * Z2 + 0.8 * R2 + 0.3 * V))
    rbind(
      data.frame(id = 1:n, time = 1, V = V, A = A1, R = R1, Z = Z1,
                 Y = NA_real_, lag_A = 0),
      data.frame(id = 1:n, time = 2, V = V, A = A2, R = R2, Z = Z2,
                 Y = Y2, lag_A = A1)
    )
  }
  # Exact truth for J^{a, a'}: R_2 under regime a (A1 = A2 = a), mediator law
  # at a' given the a-world R_2, Y at (a, Z2, R2, V).
  truth_one <- function(a, ap) {
    val <- 0
    for (V in 0:1) {
      pr2 <- expit(-0.4 + 0.9 * a + 1.4 * a + 0.3 * V)
      for (r2 in 0:1) {
        pz2 <- expit(-0.3 + 0.8 * ap + 0.7 * r2 + 0.3 * V)
        for (z2 in 0:1) {
          py <- expit(-0.5 + 0.7 * a + 0.9 * z2 + 0.8 * r2 + 0.3 * V)
          val <- val + 0.5 *
            (r2 * pr2 + (1 - r2) * (1 - pr2)) *
            (z2 * pz2 + (1 - z2) * (1 - pz2)) * py
        }
      }
    }
    val
  }
  truth <- c(Phi00 = truth_one(0, 0), Phi10 = truth_one(1, 0),
             Phi01 = truth_one(0, 1), Phi11 = truth_one(1, 1))

  d <- gen_lag_data(30000, seed = 42)
  fit <- suppressWarnings(mediation(
    data = d, id_var = "id", base_vars = "V", exposure = "A",
    outcome = "Y", time_var = "time",
    models = list(
      spec_model(A ~ V + lag_A,     var_type = "binary", mod_type = "exposure"),
      spec_model(R ~ V + A + lag_A, var_type = "binary", mod_type = "covariate"),
      spec_model(Z ~ V + A + R,     var_type = "binary", mod_type = "mediator"),
      spec_model(Y ~ V + A + Z + R, var_type = "binary", mod_type = "outcome")
    ),
    init_recode = recodes(lag_A = 0),
    in_recode   = recodes(lag_A = A),
    mediation_type = "N", estimator = "tmle", R = 1, quiet = TRUE
  ))

  est <- setNames(fit$effect_size$Est, fit$effect_size$Intervention)
  sds <- setNames(fit$effect_size$Sd,  fit$effect_size$Intervention)
  for (ph in names(truth)) {
    testthat::expect_lt(abs(est[[ph]] - truth[[ph]]), 4 * sds[[ph]] + 0.01,
                        label = ph)
  }
})


testthat::test_that("tmle cross-checks the medltmle survival + censoring + intermediate-confounder DGP", {
  testthat::skip_on_cran()

  # Compact single-fit gate on the DGP from medltmle (Malenica, Zheng &
  # van der Laan; github.com/imalenica/medltmle), the reference longitudinal-
  # mediation TMLE by the paper's own authors. This is the only test that
  # exercises the survival + censoring + intermediate-confounder combination
  # together. Large-sample truths (10^6 intervened draws, seeded) live in
  # references/validation/medltmle-crosscheck.R; a fuller bias/coverage study
  # is saved there. Here: one n = 1000 fit, checked for sane, on-target output.
  expit <- stats::plogis
  gen <- function(n, seed) {
    set.seed(seed)
    W1 <- rbinom(n, 1, 0.4); W2 <- rbinom(n, 1, 0.6)
    alive <- rep(TRUE, n); incare <- rep(TRUE, n)
    A_prev <- rep(0, n); LA_prev <- rep(0, n); LZ_prev <- rep(0, n)
    rows <- vector("list", 2L)
    for (t in 1:2) {
      lagA  <- if (t == 1) rep(0, n) else A_prev
      lagLA <- if (t == 1) rep(0, n) else LA_prev
      lagLZ <- if (t == 1) rep(0, n) else LZ_prev
      atrisk <- alive & incare
      # medltmle's C is P(uncensored); causalMed uses C = 1 => censored.
      punc <- if (t == 1) expit(3 * W2 - 0.2 * W1)
              else        expit(3 * W2 - 0.5 * lagA - 1.4 * lagLZ)
      unc  <- rbinom(n, 1, punc)
      pa   <- if (t == 1) expit(-0.1 + 0.7 * W1 + 1.2 * W2)
              else        expit(-0.1 + 1.2 * W2 + 0.7 * lagLZ - 0.1 * lagA)
      A    <- rbinom(n, 1, pa); incare <- incare & (unc == 1)
      pla  <- if (t == 1) expit(-0.8 + 0.1 * W1 + 0.3 * W2 + A)
              else        expit(-0.8 + 0.1 * lagLZ + 0.3 * lagLA + A)
      LA   <- rbinom(n, 1, pla)
      Z    <- rbinom(n, 1, expit(-0.5 + 0.8 * W2 + 0.8 * A + 1 * LA))
      plz  <- if (t == 1) expit(-1 + 0.3 * W2 + A + 0.7 * Z)
              else        expit(-1 + 0.3 * W2 + A + 0.7 * Z - 0.2 * lagLZ)
      LZ   <- rbinom(n, 1, plz)
      py   <- if (t == 1) expit(-6 + 1.5 * W2 + LA + 0.2 * LZ - 0.3 * A - 0.3 * Z - 0.2 * A * Z)
              else        expit(-6 + 1.5 * W2 + LA + 0.2 * LZ - 0.3 * A - 0.3 * Z - 0.2 * A * Z - 0.1 * lagLA)
      Y    <- rbinom(n, 1, py)
      LA[!incare] <- 0; Z[!incare] <- 0; LZ[!incare] <- 0; Y[!incare] <- 0; A[!incare] <- 0
      rows[[t]] <- data.frame(id = seq_len(n), time = t, W1 = W1, W2 = W2,
                              C = 1L - unc, A = A, LA = LA, Z = Z, LZ = LZ, Y = Y,
                              lag_A = lagA, lag_LA = lagLA, lag_LZ = lagLZ, atrisk = atrisk)
      A_prev <- A; LA_prev <- LA; LZ_prev <- LZ
      alive  <- alive & (Y == 0 | !incare)
    }
    out <- rbind(rows[[1]][rows[[1]]$atrisk, ], rows[[2]][rows[[2]]$atrisk, ])
    out$atrisk <- NULL
    out[order(out$id, out$time), ]
  }
  # Large-sample truth (seeded 10^6 draws; medltmle-crosscheck.R).
  truth <- c(Phi00 = 0.0213, Phi10 = 0.0181, Phi01 = 0.0204, Phi11 = 0.0172)

  models <- list(
    spec_model(C  ~ W2 + I((time == 1) * W1) + lag_A + lag_LZ,
               var_type = "binary", mod_type = "censor"),
    spec_model(A  ~ W2 + I((time == 1) * W1) + lag_LZ + lag_A,
               var_type = "binary", mod_type = "exposure"),
    spec_model(LA ~ A + I((time == 1) * W1) + I((time == 1) * W2) + lag_LZ + lag_LA,
               var_type = "binary", mod_type = "covariate"),
    spec_model(Z  ~ W2 + A + LA,           var_type = "binary", mod_type = "mediator"),
    spec_model(LZ ~ W2 + A + Z + lag_LZ,   var_type = "binary", mod_type = "covariate"),
    spec_model(Y  ~ W2 + LA + LZ + A + Z + A:Z + lag_LA,
               var_type = "binary", mod_type = "survival")
  )

  fit <- suppressWarnings(mediation(
    data = gen(1000, seed = 1), id_var = "id", base_vars = c("W1", "W2"),
    exposure = "A", outcome = "Y", time_var = "time", models = models,
    init_recode = recodes(lag_A = 0, lag_LA = 0, lag_LZ = 0),
    in_recode   = recodes(lag_A = A, lag_LA = LA, lag_LZ = LZ),
    mediation_type = "N", estimator = "tmle", R = 1, quiet = TRUE
  ))

  est <- setNames(fit$effect_size$Est, fit$effect_size$Intervention)
  sds <- setNames(fit$effect_size$Sd,  fit$effect_size$Intervention)
  # Sane risks and on-target functionals (self-calibrated tolerance; the fit is
  # noisy at n = 1000 with this rare outcome, so allow 4 IC-SE plus slack).
  testthat::expect_true(all(is.finite(est)) && all(est > 0 & est < 1))
  for (ph in names(truth))
    testthat::expect_lt(abs(est[[ph]] - truth[[ph]]), 4 * sds[[ph]] + 0.02, label = ph)

  # Exact natural-effect decomposition NDE + NIE = TE.
  rd <- setNames(fit$estimate$RD, fit$estimate$Effect)
  testthat::expect_equal(rd[["Direct effect"]] + rd[["Indirect effect"]],
                         rd[["Total effect"]], tolerance = 1e-10)

  # LA is exposure-affected -> the intermediate-confounder path must be flagged.
  testthat::expect_true("LA" %in% fit$intermediate_confounders)
})
