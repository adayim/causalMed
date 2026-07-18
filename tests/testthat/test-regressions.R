
testthat::test_that("spec_model recode is applied before model fitting", {
  data("nonsurvivaldata", package = "causalMed")
  dat <- data.table::as.data.table(nonsurvivaldata)
  # Column exists but holds the wrong values; the model recode must overwrite
  # it with V^2 BEFORE the outcome model is fitted.
  dat[, V_sq := 0]

  mod1 <- spec_model(A ~ V, var_type = "binary", mod_type = "exposure")
  mod3 <- spec_model(Y_bin ~ A + V_sq, var_type = "binary", mod_type = "outcome",
                     recode = recodes(V_sq = V^2))

  fit <- suppressWarnings(gformula(
    data = dat, id_var = "id", base_vars = "V",
    exposure = "A", time_var = "time",
    models = list(mod1, mod3),
    intervention = list(always = 1),
    mc_sample = 200, R = 0, quiet = TRUE, seed = 1,
    return_fitted = TRUE
  ))
  got <- coef(fit$fitted_models[["Y_bin"]])[["V_sq"]]

  # Reference: fit the same model with the recode applied by hand.
  dat2 <- data.table::copy(dat)[, V_sq := V^2]
  want <- coef(glm(Y_bin ~ A + V_sq, family = binomial(), data = dat2))[["V_sq"]]

  testthat::expect_false(is.na(got))
  testthat::expect_equal(got, want, tolerance = 1e-8)
})


testthat::test_that("per-model recode does not leak into the caller's data", {
  data("nonsurvivaldata", package = "causalMed")
  dat <- data.table::as.data.table(nonsurvivaldata)
  dat[, V_sq := 0]

  mod1 <- spec_model(A ~ V, var_type = "binary", mod_type = "exposure")
  mod3 <- spec_model(Y_bin ~ A + V_sq, var_type = "binary", mod_type = "outcome",
                     recode = recodes(V_sq = V^2))

  suppressWarnings(gformula(
    data = dat, id_var = "id", base_vars = "V",
    exposure = "A", time_var = "time",
    models = list(mod1, mod3),
    intervention = list(always = 1),
    mc_sample = 200, R = 0, quiet = TRUE, seed = 1
  ))
  # The user's data.table must not have been modified in place.
  testthat::expect_true(all(dat$V_sq == 0))
})


testthat::test_that("logical static interventions behave like numeric {0,1}", {
  data("nonsurvivaldata", package = "causalMed")

  mod1 <- spec_model(A ~ V, var_type = "binary", mod_type = "exposure")
  mod3 <- spec_model(Y_bin ~ A + V, var_type = "binary", mod_type = "outcome")
  models <- list(mod1, mod3)

  fit_log <- suppressWarnings(gformula(
    data = nonsurvivaldata, id_var = "id", base_vars = "V",
    exposure = "A", time_var = "time", models = models,
    intervention = list(natural = NULL, always = TRUE),
    mc_sample = 500, R = 0, quiet = TRUE, seed = 1
  ))
  fit_num <- suppressWarnings(gformula(
    data = nonsurvivaldata, id_var = "id", base_vars = "V",
    exposure = "A", time_var = "time", models = models,
    intervention = list(natural = NULL, always = 1),
    mc_sample = 500, R = 0, quiet = TRUE, seed = 1
  ))
  testthat::expect_equal(fit_log$effect_size$Est, fit_num$effect_size$Est)
})


testthat::test_that("more than one NULL intervention gives a clear error", {
  data("nonsurvivaldata", package = "causalMed")

  mod1 <- spec_model(A ~ V, var_type = "binary", mod_type = "exposure")
  mod3 <- spec_model(Y_bin ~ A + V, var_type = "binary", mod_type = "outcome")

  testthat::expect_error(
    gformula(
      data = nonsurvivaldata, id_var = "id", base_vars = "V",
      exposure = "A", time_var = "time", models = list(mod1, mod3),
      intervention = list(nat1 = NULL, nat2 = NULL),
      mc_sample = 100, R = 0, quiet = TRUE
    ),
    regexp = "Only one intervention element may be NULL"
  )
})


testthat::test_that("mediation(quiet = TRUE) disables the bootstrap progress bar", {
  data("nonsurvivaldata", package = "causalMed")

  captured <- NULL
  fake_arms <- list(nat0 = 0.10, nat1 = 0.20,
                    Phi00 = 0.11, Phi10 = 0.15, Phi11 = 0.19)
  testthat::local_mocked_bindings(
    bootstrap_helper = function(..., progress_bar = TRUE) {
      captured <<- progress_bar
      list(list(gform.data = fake_arms), list(gform.data = fake_arms))
    },
    .package = "causalMed"
  )

  m_med <- spec_model(M     ~ V + A + time, var_type = "normal", mod_type = "mediator")
  m_Y   <- spec_model(Y_bin ~ V + A + M,    var_type = "binary", mod_type = "outcome")

  suppressWarnings(mediation(
    data = nonsurvivaldata, id_var = "id", base_vars = "V",
    exposure = "A", outcome = "Y_bin", time_var = "time",
    models = list(m_med, m_Y),
    mediation_type = "I",
    mc_sample = 200, R = 2, quiet = TRUE, seed = 1
  ))

  testthat::expect_false(captured)
})


testthat::test_that("mediation recovers the documented Yamamuro true values on yamamurodata", {
  testthat::skip_on_cran()
  data("yamamurodata", package = "causalMed")

  models_yam <- list(
    spec_model(A ~ V + lag1_A + lag1_L + I(lag1_L^2) + lag1_M1 + lag1_M2 +
                 I(as.integer(time == 2)),
               var_type = "binary", mod_type = "exposure", subset = time > 0),
    spec_model(L ~ V + A + lag1_L + I(lag1_L^2) + lag1_M1 + lag1_M2 +
                 I(as.integer(time == 2)),
               var_type = "normal", mod_type = "covariate", subset = time > 0),
    spec_model(M1 ~ V + A + L + I(L^2) + lag1_L + I(lag1_L^2) + lag1_M1 + lag1_M2 +
                 I(as.integer(time == 2)),
               var_type = "normal", mod_type = "mediator", subset = time > 0),
    spec_model(M2 ~ V + A + L + I(L^2) + lag1_L + I(lag1_L^2) + lag1_M1 + lag1_M2 +
                 I(as.integer(time == 2)),
               var_type = "normal", mod_type = "mediator", subset = time > 0),
    spec_model(Y ~ V + A + L + I(L^2) + M1 + M2 + M1:M2 +
                 I(as.integer(time == 1)) + I(as.integer(time == 2)),
               var_type = "binary", mod_type = "survival")
  )

  fit <- suppressWarnings(mediation(
    data = yamamurodata, id_var = "id", time_var = "time",
    base_vars = c("V", "L0base", "M10base", "M20base"),
    exposure = "A", outcome = "Y", models = models_yam,
    init_recode = recodes(L = L0base, M1 = M10base, M2 = M20base,
                          lag1_A = 0, lag1_L = 0, lag1_M1 = 0, lag1_M2 = 0),
    in_recode = recodes(lag1_A = A, lag1_L = L, lag1_M1 = M1, lag1_M2 = M2),
    mediation_type = "I", n_vw = 1,
    mc_sample = 5000, R = 1, quiet = TRUE, seed = 7
  ))

  est <- fit$estimate
  pick <- function(lbl) est$RD[est$Effect == lbl] * 100

  # Documented large-sample true values (?yamamurodata), in percentage points.
  # Tolerance covers single-dataset sampling error plus MC noise; this guards
  # against gross regressions (sign flips, intervention mix-ups), not fine accuracy.
  truth <- c("Total effect" = -6.36, "Direct effect" = -3.20,
             "Indirect effect (M1)" = -2.29, "Indirect effect (M2)" = -0.97,
             "TE - (Direct + Indirect)" = 0.10)
  for (lbl in names(truth)) {
    testthat::expect_lt(abs(pick(lbl) - truth[[lbl]]), 2.5, label = lbl)
  }

  # Exact decomposition identity: IDE + sum(IIE) + residual == TE
  testthat::expect_equal(
    pick("Direct effect") + pick("Indirect effect (M1)") +
      pick("Indirect effect (M2)") + pick("TE - (Direct + Indirect)"),
    pick("Total effect"),
    tolerance = 1e-8
  )
})


testthat::test_that("observed_benchmark matches survival::survfit KM under censoring (gvhd)", {
  testthat::skip_if_not_installed("survival")

  data("gvhd", package = "causalMed")

  ob <- observed_benchmark(gvhd, "d", "day", is_survival = TRUE)

  # Gold standard: Kaplan-Meier cumulative incidence at end of follow-up,
  # from one record per subject (last day, ever-death status). gvhd has
  # real loss-to-follow-up (censlost) plus administrative censoring, so
  # this exercises the product-limit correction on genuinely censored data.
  per_id <- data.frame(
    time   = tapply(gvhd$day, gvhd$id, max),
    status = tapply(gvhd$d,   gvhd$id, max)
  )
  km <- survival::survfit(survival::Surv(time, status) ~ 1, data = per_id)
  km_ci <- 1 - min(km$surv)

  testthat::expect_equal(ob$value, km_ci, tolerance = 1e-12)

  # And the naive ever-event proportion must differ (censoring bias),
  # guarding against a regression to the naive estimator.
  naive <- sum(per_id$status) / nrow(per_id)
  testthat::expect_gt(abs(ob$value - naive), 0.01)
})


testthat::test_that("spec_model(truncate=) controls range clipping of normal draws", {
  # Review finding 3.3: "normal" draws are clipped to the observed range by
  # default (i.e. a bounded normal, unlike gfoRmula's plain normal).
  # truncate = FALSE must restore the untruncated Gaussian draw.
  set.seed(4)
  n <- 1500
  d <- data.frame(id = rep(seq_len(n), each = 2), time = rep(0:1, n))
  d$V <- rnorm(nrow(d))
  d$L <- rnorm(nrow(d), 0.5 * d$V, 1)
  d$A <- rbinom(nrow(d), 1, plogis(0.3 * d$L))
  d$Y <- rbinom(nrow(d), 1, plogis(-1 + 0.5 * d$A + 0.3 * d$L))
  obs_range <- range(d$L)

  sim_L_range <- function(truncate) {
    fit <- suppressWarnings(gformula(
      data = d, id_var = "id", base_vars = "V", exposure = "A",
      time_var = "time",
      models = list(
        spec_model(L ~ V + time, var_type = "normal", mod_type = "covariate",
                   truncate = truncate),
        spec_model(A ~ L + V + time, var_type = "binary", mod_type = "exposure"),
        spec_model(Y ~ A + L + V, var_type = "binary", mod_type = "outcome")
      ),
      intervention = list(always = 1),
      mc_sample = 20000, R = 1, quiet = TRUE, seed = 3, return_data = TRUE
    ))
    range(fit$sim_data$L)
  }

  r_on  <- sim_L_range(TRUE)
  r_off <- sim_L_range(FALSE)

  # Default clips into the observed support ...
  testthat::expect_gte(r_on[1], obs_range[1] - 1e-9)
  testthat::expect_lte(r_on[2], obs_range[2] + 1e-9)
  # ... while truncate = FALSE is allowed to leave it.
  testthat::expect_true(r_off[1] < obs_range[1] || r_off[2] > obs_range[2])

  # Back-compat: the flag defaults to TRUE (historical behaviour).
  testthat::expect_true(
    isTRUE(spec_model(L ~ V, var_type = "normal", mod_type = "covariate")$truncate)
  )
  testthat::expect_error(
    spec_model(L ~ V, var_type = "normal", mod_type = "covariate", truncate = NA),
    "single TRUE or FALSE"
  )
})


testthat::test_that("print/summary work when mediation_type is left at its default (review 2026-07-15 finding 2)", {
  # all.args is captured before match.arg(), so without the write-back a
  # defaulted call stored c("I", "N") and print() errored on the length-2
  # condition ("the condition has length > 1").
  data("nonsurvivaldata", package = "causalMed")

  m_med <- spec_model(M     ~ V + A + time, var_type = "normal", mod_type = "mediator")
  m_Y   <- spec_model(Y_bin ~ V + A + M,    var_type = "binary", mod_type = "outcome")

  fit <- suppressWarnings(mediation(
    data = nonsurvivaldata, id_var = "id", base_vars = "V",
    exposure = "A", outcome = "Y_bin", time_var = "time",
    models = list(m_med, m_Y),
    # mediation_type deliberately omitted -> default "I"
    mc_sample = 2000, R = 1, quiet = TRUE, seed = 1
  ))

  testthat::expect_identical(fit$all.args$mediation_type, "I")
  testthat::expect_output(print(fit), "Interventional effects")
  testthat::expect_no_error(utils::capture.output(summary(fit)))
})


testthat::test_that("a numeric ref_int is resolved to its name (no self-contrast row)", {
  # Review finding G-1: risk_estimate_point() drops the reference via
  # setdiff(names(intervention), ref_int); a numeric index never matches a
  # name, so the reference used to be contrasted against itself (RD 0, RR 1).
  data("nonsurvivaldata", package = "causalMed")

  models <- list(
    spec_model(A     ~ V,         var_type = "binary", mod_type = "exposure"),
    spec_model(Y_bin ~ A + V,     var_type = "binary", mod_type = "outcome")
  )
  ints <- list(never = 0, always = 1)

  fit_num <- suppressWarnings(gformula(
    data = nonsurvivaldata, id_var = "id", base_vars = "V", exposure = "A",
    time_var = "time", models = models, intervention = ints, ref_int = 1,
    mc_sample = 2000, R = 1, quiet = TRUE, seed = 9
  ))

  # ref_int = 1 means "never": only the "always" contrast should be reported.
  testthat::expect_equal(fit_num$all.args$ref_int, "never")
  testthat::expect_false(any(grepl("never - never|never / never",
                                   fit_num$estimate$Intervention)))
  testthat::expect_setequal(fit_num$estimate$Intervention,
                            c("always - never", "always / never"))

  # ... and it must agree with passing the reference by name.
  fit_chr <- suppressWarnings(gformula(
    data = nonsurvivaldata, id_var = "id", base_vars = "V", exposure = "A",
    time_var = "time", models = models, intervention = ints, ref_int = "never",
    mc_sample = 2000, R = 1, quiet = TRUE, seed = 9
  ))
  testthat::expect_equal(fit_num$estimate$Estimate, fit_chr$estimate$Estimate)
})
