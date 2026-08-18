# Regression tests for the package.
#
# Each block names a defect found in a past review, so that a change which
# reintroduces the defect fails loudly. 

# ---- Shared fixtures --------------------------------------------------------

# spec_model() captures its call unevaluated, so the formula must be written
# out literally -- passing it through a variable would store the symbol.
make_models <- function() {
  list(
    spec_model(A ~ V + lag1_A + time, var_type = "binary", mod_type = "exposure"),
    spec_model(L1 ~ V + A + lag1_L1 + time, var_type = "normal", mod_type = "covariate"),
    spec_model(Y_bin ~ V + A + L1, var_type = "binary", mod_type = "outcome")
  )
}

# Same, but with `time` in the outcome model -- an aliased term, because Y_bin
# is recorded only at the last time point.
make_models_rank_deficient <- function() {
  list(
    spec_model(A ~ V + lag1_A + time, var_type = "binary", mod_type = "exposure"),
    spec_model(L1 ~ V + A + lag1_L1 + time, var_type = "normal", mod_type = "covariate"),
    spec_model(Y_bin ~ V + A + L1 + time, var_type = "binary", mod_type = "outcome")
  )
}

make_models_no_exposure <- function() {
  list(
    spec_model(L1 ~ V + A + lag1_L1 + time, var_type = "normal", mod_type = "covariate"),
    spec_model(Y_bin ~ V + A + L1, var_type = "binary", mod_type = "outcome")
  )
}

run_gf <- function(models, ...) {
  gformula(data = nonsurvivaldata, id_var = "id", base_vars = "V",
           exposure = "A", time_var = "time", models = models,
           init_recode = recodes(lag1_A = 0, lag1_L1 = 0),
           in_recode   = recodes(lag1_A = A, lag1_L1 = L1),
           mc_sample = 400, R = 1, quiet = TRUE, seed = 1, ...)
}

# Shared DGP for the two interventional mediator-pool tests. The mediator
# climbs by ~10 units per period (means near 0, 10, 20, 30), so any collapse or
# shift of the per-time pool slices is unmistakable rather than subtle.
make_pool_data <- function(tt = 0:3, n = 600) {
  set.seed(4242)
  d <- data.table::CJ(id = seq_len(n), time = tt)
  data.table::setorderv(d, c("id", "time"))
  d[, V := rep(stats::rnorm(n), each = length(tt))]
  d[, A := stats::rbinom(.N, 1, 0.5)]
  d[, L := V + 0.5 * A + stats::rnorm(.N, 0, 0.5)]
  d[, M := 10 * time + 2 * A + L + stats::rnorm(.N, 0, 1)]
  d[, Y := stats::rbinom(.N, 1, stats::plogis(-2 + 0.5 * A + 0.02 * M + 0.3 * L))]
  d[, lag1_A := 0]; d[, lag1_M := 0]; d[, lag1_L := 0]
  d[]
}

# Monte Carlo cohort for the pool tests: id + baseline only, as the engine
# builds it. Called after make_pool_data() so the RNG stream is deterministic.
make_pool_mc <- function(d, size = 800) {
  base_dat <- unique(d[, c("id", "V"), with = FALSE])
  data.table::as.data.table(base_dat[sample(nrow(base_dat), size, replace = TRUE), ])
}


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
  # risk_estimate_point() drops the reference via
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

testthat::test_that("ref_int = 'natural' creates the natural course when it was not supplied", {
  data(nonsurvivaldata)

  # Used to sail past every check and yield a zero-row contrast table.
  fit <- run_gf(make_models(), intervention = list(always = 1, never = 0),
                ref_int = "natural")

  testthat::expect_true("natural" %in% as.character(fit$effect_size$Intervention))
  testthat::expect_equal(nrow(fit$estimate), 4L)   # 2 contrasts x {RD, RR}
  testthat::expect_true(all(is.finite(fit$estimate$Estimate)))
})


testthat::test_that("ref_int = 'natural' reuses a NULL element supplied under another name", {
  data(nonsurvivaldata)

  # `ref_int = "natural"` means "contrast against the natural course", and the
  # natural course is the NULL element -- whatever the user called it. Testing
  # only for the NAME "natural" prepended a SECOND NULL element here, breaking
  # the one-NULL-element rule that check_intervention() enforces, simulating
  # the natural course twice and emitting a nonsense `nat - natural` row.
  fit <- run_gf(make_models(), intervention = list(nat = NULL, always = 1),
                ref_int = "natural")

  interventions <- as.character(fit$effect_size$Intervention)
  testthat::expect_equal(sort(interventions), c("always", "nat"))
  testthat::expect_false("natural" %in% interventions)
  testthat::expect_equal(fit$all.args$ref_int, "nat")
  testthat::expect_equal(nrow(fit$estimate), 2L)   # 1 contrast x {RD, RR}

  # A NON-NULL intervention the user chose to call "natural" is still honoured
  # as the reference, and must not collide with an auto-added one.
  fit2 <- run_gf(make_models(), intervention = list(natural = 1, never = 0),
                 ref_int = "natural")
  testthat::expect_equal(sum(names(fit2$all.args$intervention) == "natural"), 1L)
  testthat::expect_equal(fit2$all.args$ref_int, "natural")
})


testthat::test_that("intervention = NULL is labelled 'natural'", {
  data(nonsurvivaldata)
  fit <- run_gf(make_models(), intervention = NULL)
  testthat::expect_equal(as.character(fit$effect_size$Intervention), "natural")
  testthat::expect_equal(fit$all.args$ref_int, "natural")
})


testthat::test_that("ref_int validation rejects out-of-range and non-whole numbers", {
  data(nonsurvivaldata)
  testthat::expect_error(
    run_gf(make_models(), intervention = list(always = 1), ref_int = -1),
    "whole number")
  testthat::expect_error(
    run_gf(make_models(), intervention = list(always = 1), ref_int = 1.5),
    "whole number")
  testthat::expect_error(
    run_gf(make_models(), intervention = list(always = 1), ref_int = 99),
    "must be less than the length of intervention")
})


testthat::test_that("a natural course without an exposure model fails with an actionable message", {
  data(nonsurvivaldata)

  testthat::expect_error(run_gf(make_models_no_exposure(), intervention = NULL),
                         "no exposure model was supplied")

  # With the default ref_int = 0 the natural course is auto-prepended as the
  # reference, so a "static-only" list still needs an exposure model. The
  # message must therefore NOT advise "use only static/dynamic interventions"
  # -- that is a no-op which lands the user back on this same error. The
  # actionable remedy is to name an existing intervention as `ref_int`.
  err <- testthat::expect_error(
    run_gf(make_models_no_exposure(), intervention = list(always = 1, never = 0)),
    "no exposure model was supplied")
  testthat::expect_match(conditionMessage(err), "ref_int")

  # Naming a static intervention as the reference means no natural course is
  # added, and then no exposure model is required.
  testthat::expect_s3_class(
    run_gf(make_models_no_exposure(), intervention = list(always = 1, never = 0),
           ref_int = "never"),
    "gformula")
})


testthat::test_that("a rank-deficient model drops aliased terms instead of returning all-NA", {
  data(nonsurvivaldata)

  # `time` is constant among the rows where Y_bin is observed (end of
  # follow-up only), so glm() aliases it. The cached `beta` then held an NA
  # and `mm %*% beta` made EVERY simulated value NA.
  msgs <- character(0)
  fit <- withCallingHandlers(
    run_gf(make_models_rank_deficient(),
           intervention = list(natural = NULL, always = 1)),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
  )

  testthat::expect_true(all(is.finite(fit$effect_size$Est)))
  testthat::expect_true(any(grepl("rank deficient", msgs)))
  testthat::expect_true(any(grepl("time", msgs[grepl("rank deficient", msgs)])))

  # The retained coefficients must reproduce predict()'s linear predictor.
  g  <- stats::glm(Y_bin ~ V + A + L1 + time, family = stats::binomial(),
                   data = nonsurvivaldata)
  nd <- data.frame(V = 0.5, A = 1, L1 = 0.2, time = 3)
  b  <- stats::coef(g); b <- b[!is.na(b)]
  mm <- stats::model.matrix(stats::delete.response(stats::terms(g)), nd)
  testthat::expect_equal(
    as.numeric(stats::plogis(mm[, names(b), drop = FALSE] %*% b)),
    suppressWarnings(as.numeric(stats::predict(g, nd, type = "response")))
  )
})


testthat::test_that("aliased coefficients are dropped by position when coef() is unnamed", {
  data(nonsurvivaldata)

  # The alias drop indexed `beta` by NAME. A fit whose coef() vector is unnamed
  # -- legal, and what a hand-rolled or wrapped learner may return -- reduced
  # `beta` to numeric(0), reported an empty `{}` dropped-terms list, and then
  # refused the model with "provides no usable model terms/coefficients".
  fit_unnamed <- function(formula, data, ...) {
    g <- stats::glm(formula, data = data, ...)
    # A rank-deficient fit whose coefficients arrive without names.
    b <- stats::coef(g)
    b[2L] <- NA_real_
    names(b) <- NULL
    g$coefficients <- b
    g
  }

  assign("causalMedTestFitUnnamed", fit_unnamed, envir = globalenv())
  on.exit(suppressWarnings(rm("causalMedTestFitUnnamed", envir = globalenv())),
          add = TRUE)

  # var_type = "custom" with no custom_sim: the simulation still needs this
  # node's linear predictor, so the alias drop has to leave `beta` usable.
  models <- list(
    spec_model(A ~ V + lag1_A + time, var_type = "binary", mod_type = "exposure"),
    spec_model(L1 ~ V + A + lag1_L1 + time, var_type = "custom",
               mod_type = "covariate", custom_fit = causalMedTestFitUnnamed),
    spec_model(Y_bin ~ V + A + L1, var_type = "binary", mod_type = "outcome")
  )

  # `causalMed:::causalmed_env$warning <- NULL` is not assignable; take the
  # environment first and mutate it by reference.
  cm_env <- causalMed:::causalmed_env
  cm_env$warning <- NULL
  fitted <- causalMed:::fit_spec_models(models, nonsurvivaldata)
  warns  <- cm_env$warning

  # The estimable coefficients survive: intercept + 4 terms, 1 aliased, 4 kept.
  testthat::expect_equal(length(fitted[[2L]]$beta), 4L)
  testthat::expect_false(anyNA(fitted[[2L]]$beta))
  # And the design matrix is subset to match, by position.
  testthat::expect_equal(fitted[[2L]]$beta_cols, c(1L, 3L, 4L, 5L))

  # The warning must name the dropped term, not report an empty `{}` set.
  rank_msg <- warns[grepl("rank deficient", warns)]
  testthat::expect_length(rank_msg, 1L)
  testthat::expect_false(grepl("\\{\\}", rank_msg))

  # lin_pred() must then produce finite values, not fail on a zero-length beta.
  lp <- causalMed:::lin_pred(
    fitted[[2L]],
    data.table::data.table(V = 0.5, A = 1, lag1_L1 = 0.2, time = 3))
  testthat::expect_true(is.finite(lp))
})


testthat::test_that("the interventional mediator pool preserves the per-time trajectory", {
  # The pool stores each reference-cohort individual's WHOLE M(1:T) path, and a
  # single permutation hands subject i one pool individual's entire path. A
  # first attempt at type-preserving pool storage kept the data.table's live
  # column vector instead of a snapshot; because simulate_data() writes the next
  # time step with `data[cond, (col) := vals]` -- a `:=` with an `i` argument,
  # which sub-assigns IN PLACE -- every stored slice mutated together and the
  # pool collapsed to M(T) repeated. That silently distorted every
  # interventional Phi .
  #
  # The pool is inspected directly -- `sim_data` cannot be used for this,
  # because it is the end-of-follow-up snapshot, not a time series.
  tt <- 0:3
  d  <- make_pool_data(tt)

  models <- list(
    spec_model(A ~ V + lag1_A + time, var_type = "binary", mod_type = "exposure"),
    spec_model(L ~ V + A + lag1_L + time, var_type = "normal", mod_type = "covariate"),
    spec_model(M ~ V + A + L + lag1_M + time, var_type = "normal", mod_type = "mediator"),
    spec_model(Y ~ V + A + M + L, var_type = "binary", mod_type = "outcome")
  )

  fit_mods <- causalMed:::fit_spec_models(models, d)
  df_mc    <- make_pool_mc(d)

  run <- causalMed:::simulate_intervention(
    data           = df_mc,
    models         = fit_mods,
    exposure       = "A",
    time_var       = "time",
    time_seq       = tt,
    intervention   = causalMed:::intervention_spec(0, list()),
    init_recode    = recodes(lag1_A = 0, lag1_M = 0, lag1_L = 0),
    in_recode      = recodes(lag1_A = A, lag1_M = M, lag1_L = L),
    mediation_type = "I",
    collect_pool   = TRUE
  )
  pool <- run$pool[["M"]]

  # Shape: one slice per time point, one entry per pool individual.
  testthat::expect_equal(length(pool), length(tt))
  testthat::expect_true(all(vapply(pool, length, integer(1)) == nrow(df_mc)))

  # Trajectory: the per-time means must track the DGP's ~10-per-period climb.
  # Under the collapsed-pool bug every slice equalled the last one, so this
  # range was 0.
  slice_means <- vapply(pool, mean, numeric(1))
  testthat::expect_gt(diff(range(slice_means)), 20)
  testthat::expect_true(all(diff(slice_means) > 5))

  # And it must be a per-individual trajectory, not one value repeated: the
  # first and last slices should differ for essentially every pool individual.
  testthat::expect_gt(mean(pool[[1L]] != pool[[length(pool)]]), 0.99)
})


testthat::test_that("the mediator pool keeps one slice per time point when a subset skips a step", {
  # simulate_data() skips a model entirely when its `subset` selects no rows
  # (`if (sum(cond) != 0L)`), so at such a time step the mediator column does
  # not exist yet and `data[[mv]]` is NULL. `m_pool_out[[mv]][[indx]] <- NULL`
  # DELETES the list element instead of storing it: the pool silently shrinks,
  # every later time slice shifts down one position, and the run finally dies
  # in sim_value() with "Argument eta must be a nonempty numeric vector".
  tt <- 0:3
  d  <- make_pool_data(tt)

  # The mediator is only defined from t = 1 on, so time 0 selects no rows and
  # the M column is never created at that step. The outcome is restricted the
  # same way -- nothing may reference M before it exists, or the run dies of
  # that instead of the pool defect under test.
  models <- list(
    spec_model(A ~ V + lag1_A + time, var_type = "binary", mod_type = "exposure"),
    spec_model(L ~ V + A + lag1_L + time, var_type = "normal", mod_type = "covariate"),
    spec_model(M ~ V + A + L + time, var_type = "normal",
               mod_type = "mediator", subset = time > 0),
    spec_model(Y ~ V + A + M + L, var_type = "binary", mod_type = "outcome",
               subset = time > 0)
  )

  fit_mods <- causalMed:::fit_spec_models(models, d)
  df_mc    <- make_pool_mc(d)

  run <- causalMed:::simulate_intervention(
    data           = df_mc,
    models         = fit_mods,
    exposure       = "A",
    time_var       = "time",
    time_seq       = tt,
    intervention   = causalMed:::intervention_spec(0, list()),
    init_recode    = recodes(lag1_A = 0, lag1_L = 0),
    in_recode      = recodes(lag1_A = A, lag1_L = L),
    mediation_type = "I",
    collect_pool   = TRUE
  )
  pool <- run$pool[["M"]]

  # One slice per time point, every slice the full pool -- no shrinking, no
  # shifting. The skipped step carries an all-missing slice, not a deletion.
  testthat::expect_equal(length(pool), length(tt))
  testthat::expect_true(all(vapply(pool, length, integer(1)) == nrow(df_mc)))
  testthat::expect_true(all(is.na(pool[[1L]])))
  testthat::expect_false(anyNA(pool[[length(tt)]]))

  # The surviving slices must still track the DGP's ~10-per-period climb, i.e.
  # slice k must really be time_seq[k] and not a shifted neighbour.
  later_means <- vapply(pool[-1L], mean, numeric(1))
  testthat::expect_true(all(diff(later_means) > 5))
})


testthat::test_that("a categorical mediator keeps its levels through the interventional pool", {
  testthat::skip_if_not_installed("Hmisc")
  data(nonsurvivaldata)

  dat <- data.table::as.data.table(data.table::copy(nonsurvivaldata))
  data.table::set(dat, j = "Mcat",
                  value = cut(dat$M, breaks = 3, labels = c("lo", "mid", "hi")))

  models <- list(
    spec_model(A ~ V + lag1_A + time, var_type = "binary", mod_type = "exposure"),
    spec_model(L1 ~ V + A + lag1_L1 + time, var_type = "normal", mod_type = "covariate"),
    spec_model(Mcat ~ V + A + L1 + time, var_type = "categorical", mod_type = "mediator"),
    spec_model(Y_bin ~ V + A + Mcat + L1, var_type = "binary", mod_type = "outcome")
  )

  # The pool used to be a numeric matrix, so a factor mediator was silently
  # reduced to its integer level codes before being fed back to the models.
  fit <- suppressWarnings(suppressMessages(
    mediation(data = dat, id_var = "id", base_vars = "V", exposure = "A",
              outcome = "Y_bin", time_var = "time", models = models,
              init_recode = recodes(lag1_A = 0, lag1_L1 = 0),
              in_recode   = recodes(lag1_A = A, lag1_L1 = L1),
              mediation_type = "I", mc_sample = 300, R = 1,
              quiet = TRUE, seed = 5, return_data = TRUE)
  ))

  testthat::expect_true(all(is.finite(fit$effect_size$Est)))
  testthat::expect_true(all(is.finite(fit$estimate$RD)))
  testthat::expect_true(
    all(as.character(unique(fit$sim_data$Mcat)) %in% c("lo", "mid", "hi"))
  )
})


testthat::test_that("base_vars that are not time-fixed are reported", {
  data(nonsurvivaldata)
  testthat::expect_warning(
    gformula(data = nonsurvivaldata, id_var = "id", base_vars = c("V", "L2"),
             exposure = "A", time_var = "time", models = make_models(),
             intervention = list(always = 1),
             init_recode = recodes(lag1_A = 0, lag1_L1 = 0),
             in_recode   = recodes(lag1_A = A, lag1_L1 = L1),
             mc_sample = 200, R = 1, quiet = TRUE, seed = 1),
    "not time-fixed"
  )
})


testthat::test_that("mediation() requires `outcome` to match the outcome model's response", {
  data(nonsurvivaldata)
  med_models <- list(
    spec_model(A ~ V + lag1_A + time, var_type = "binary", mod_type = "exposure"),
    spec_model(L1 ~ V + A + lag1_L1 + time, var_type = "normal", mod_type = "covariate"),
    spec_model(M ~ V + A + L1 + lag1_M + time, var_type = "normal", mod_type = "mediator"),
    spec_model(Y_bin ~ V + A + M + L1, var_type = "binary", mod_type = "outcome")
  )
  run_med <- function(outcome, ...) {
    mediation(data = nonsurvivaldata, id_var = "id", base_vars = "V",
              exposure = "A", outcome = outcome, time_var = "time",
              models = med_models,
              init_recode = recodes(lag1_A = 0, lag1_L1 = 0, lag1_M = 0),
              in_recode   = recodes(lag1_A = A, lag1_L1 = L1, lag1_M = M),
              mediation_type = "I", mc_sample = 300, R = 1,
              quiet = TRUE, seed = 3, ...)
  }
  # A typo used to give a silent NA benchmark, and made the TMLE binary-outcome
  # guard pass vacuously on a zero-length vector.
  testthat::expect_error(run_med("TYPO_NOT_A_COLUMN"), "must name the response")
  testthat::expect_s3_class(run_med("Y_bin"), "gformula")
})


testthat::test_that("spec_model validates custom_sim and flags an unused custom_fit", {
  testthat::expect_error(
    spec_model(M ~ A, var_type = "custom", custom_sim = "not a function"),
    "must be a function")
  testthat::expect_error(
    spec_model(M ~ A, var_type = "custom", custom_sim = function(fit) fit),
    "two arguments")
  testthat::expect_warning(
    spec_model(M ~ A, var_type = "normal", custom_fit = stats::glm),
    "only used when")
  # A well-formed custom_sim is accepted.
  testthat::expect_s3_class(
    spec_model(M ~ A, var_type = "custom",
               custom_sim = function(fit, newdf) rep(0, nrow(newdf))),
    "causalMed_gmodel")
})


testthat::test_that("non-spec_model entries are rejected before any field is read", {
  data(nonsurvivaldata)
  testthat::expect_error(
    gformula(data = nonsurvivaldata, id_var = "id", base_vars = "V",
             exposure = "A", time_var = "time",
             models = list(spec_model(Y_bin ~ A, var_type = "binary",
                                      mod_type = "outcome"),
                           list(mod_type = "covariate"))),
    "must be `causalMed_gmodel` object"
  )
})


testthat::test_that("a custom_fit model of an arbitrary class works behind custom_sim", {
  data(nonsurvivaldata)

  # Stands in for a data-adaptive learner (ranger/xgboost/SuperLearner): no
  # `terms` component, no coef() method, its own predict(). Previously
  # fit_spec_models() died extracting terms()/coef(), and even when a fit did
  # survive extraction, wrap_fitted_models() then crashed on
  # summary(fit)$coefficients -- but only with the DEFAULT return_fitted =
  # FALSE, so a valid model failed on the default arguments.
  # spec_model() stores `custom_fit` as an unevaluated SYMBOL and the fit is
  # evaluated inside the package namespace, so the fitting function has to be
  # reachable from the global environment (or be package-qualified) -- the same
  # constraint the docs state for parallel bootstrap. A function local to this
  # test block would not be found, so publish it for the duration of the test.
  fit_bb <- function(formula, data, ...) {
    mf <- stats::model.frame(formula, data)
    y  <- stats::model.response(mf)
    structure(list(ybar = mean(y, na.rm = TRUE), sd = stats::sd(y, na.rm = TRUE)),
              class = "causalMedTestBlackbox")
  }
  registerS3method("predict", "causalMedTestBlackbox",
                   function(object, newdata, ...) rep(object$ybar, nrow(newdata)))
  # custom_sim must DRAW, not return the conditional mean.
  sim_bb <- function(fit, newdf) stats::rnorm(nrow(newdf), predict(fit, newdf), fit$sd)

  with_global_fitter <- function(code) {
    assign("causalMedTestFitBB", fit_bb, envir = globalenv())
    on.exit(suppressWarnings(rm("causalMedTestFitBB", envir = globalenv())),
            add = TRUE)
    force(code)
  }

  with_global_fitter({
    models <- list(
      spec_model(A ~ V + lag1_A + time, var_type = "binary", mod_type = "exposure"),
      spec_model(L1 ~ V + A + lag1_L1 + time, var_type = "custom",
                 mod_type = "covariate",
                 custom_fit = causalMedTestFitBB, custom_sim = sim_bb),
      spec_model(Y_bin ~ V + A + L1, var_type = "binary", mod_type = "outcome")
    )
    fit <- run_gf(models, intervention = list(natural = NULL, always = 1))
    testthat::expect_s3_class(fit, "gformula")
    testthat::expect_true(all(is.finite(fit$effect_size$Est)))

    # print()/summary() must survive the missing coefficient table and say why.
    out <- utils::capture.output(print(fit, models = TRUE))
    testthat::expect_true(any(grepl("no coefficient table available", out)))
    testthat::expect_no_error(utils::capture.output(summary(fit)))

    # Without custom_sim the linear predictor is required, so refuse early with
    # a message that names the variable and the way out -- not the bare
    # "requires numeric/complex matrix/vector arguments" from `mm %*% NULL`.
    models_no_sim <- list(
      spec_model(A ~ V + lag1_A + time, var_type = "binary", mod_type = "exposure"),
      spec_model(L1 ~ V + A + lag1_L1 + time, var_type = "custom",
                 mod_type = "covariate", custom_fit = causalMedTestFitBB),
      spec_model(Y_bin ~ V + A + L1, var_type = "binary", mod_type = "outcome")
    )
    testthat::expect_error(
      run_gf(models_no_sim, intervention = list(natural = NULL, always = 1)),
      "custom_sim")
  })
})


testthat::test_that("survival runs with out_recode emit no data.table shallow-copy warning", {
  data(survivaldata)
  dat <- data.table::as.data.table(data.table::copy(survivaldata))
  data.table::set(dat, j = "cens", value = 0L)

  models <- list(
    spec_model(A ~ V + lag1_A + time, var_type = "binary", mod_type = "exposure"),
    spec_model(cens ~ V + time, var_type = "binary", mod_type = "censor"),
    spec_model(L ~ V + A + lag1_L + time, var_type = "normal", mod_type = "covariate"),
    spec_model(Y ~ V + A + L + time, var_type = "binary", mod_type = "survival")
  )

  warns <- character(0)
  withCallingHandlers(
    gformula(data = dat, id_var = "id", base_vars = "V", exposure = "A",
             time_var = "time", models = models,
             intervention = list(always = 1),
             init_recode = recodes(lag1_A = 0, lag1_L = 0),
             in_recode   = recodes(lag1_A = A, lag1_L = L),
             out_recode  = recodes(post_treat = A * 2),
             mc_sample = 200, R = 1, quiet = TRUE, seed = 1),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  testthat::expect_false(any(grepl("shallow copy", warns)))
})


testthat::test_that("documented dataset conventions hold", {
  data(nonsurvivaldata); data(survivaldata)
  n <- nonsurvivaldata; s <- survivaldata

  # ?nonsurvivaldata: 15,000 rows / 3000 subjects; outcomes only at t = 4;
  # lag columns NA at t = 0.
  testthat::expect_equal(nrow(n), 15000L)
  testthat::expect_equal(length(unique(n$id)), 3000L)
  testthat::expect_true(all(is.na(n$Y_bin[n$time < 4])))
  testthat::expect_false(anyNA(n$Y_bin[n$time == 4]))
  testthat::expect_true(all(is.na(n$Y_cont[n$time < 4])))
  testthat::expect_true(all(is.na(n$lag1_A[n$time == 0])))

  # ?survivaldata: at-risk format, no missing values, lags 0-filled at t = 0.
  testthat::expect_equal(nrow(s), 7113L)
  testthat::expect_equal(length(unique(s$id)), 3000L)
  testthat::expect_false(anyNA(s))
  testthat::expect_true(all(s$lag1_A[s$time == 0] == 0))
})
