
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
  # against gross regressions (sign flips, arm mix-ups), not fine accuracy.
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
