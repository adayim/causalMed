# mediation_type = "N" draws the cross-world mediator from its model evaluated
# on the exposure HISTORY set to the other regime (Zheng & van der Laan 2017,
# Eq. 5), not on the current exposure alone. The swap sets exactly the exposure
# and the first-order in_recode lags; every other exposure-derived input is
# rejected up front (the check fails closed).

# ---- check_natural_exposure_history() ---------------------------------------

testthat::test_that("check_natural_exposure_history accepts what the swap sets", {
  f   <- causalMed:::check_natural_exposure_history
  out <- spec_model(Y ~ A + M, var_type = "binary", mod_type = "outcome")
  init_rc <- recodes(lag1_A = 0, lag1_M = 0, lag1_L = 0, lag2_L = 0)
  in_rc   <- recodes(lag1_A = A, lag1_M = M)

  testthat::expect_identical(
    f(list(spec_model(M ~ A + lag1_A + lag1_M, var_type = "normal", mod_type = "mediator"), out),
      "A", init_rc, in_rc, NULL),
    "lag1_A")
  # Exposure terms written in the formula are evaluated on the swapped values.
  testthat::expect_no_error(
    f(list(spec_model(M ~ A:lag1_M + I(A * lag1_A), var_type = "normal", mod_type = "mediator"), out),
      "A", init_rc, in_rc, NULL))
  # A recode and a subset that never touch the exposure, including an
  # order-dependent chain.
  testthat::expect_no_error(
    f(list(spec_model(M ~ lag2_L + V, var_type = "normal", mod_type = "mediator",
                      recode = recodes(lag2_L = lag1_L, lag1_L = L), subset = time > 0), out),
      "A", init_rc, in_rc, NULL))
  testthat::expect_identical(
    f(list(spec_model(M ~ A, var_type = "normal", mod_type = "mediator"), out), "A"),
    character(0))
})

testthat::test_that("check_natural_exposure_history rejects every other exposure-derived input", {
  f   <- causalMed:::check_natural_exposure_history
  out <- spec_model(Y ~ A + M, var_type = "binary", mod_type = "outcome")
  init_rc <- recodes(lag1_A = 0, lag2_A = 0, cumA = 0)
  in_rc   <- recodes(lag2_A = lag1_A, lag1_A = A, cumA = cumA + A)
  rejects <- function(models, pattern, init = init_rc, inr = in_rc, outr = NULL)
    testthat::expect_error(f(models, "A", init, inr, outr), pattern, fixed = TRUE)

  # chained lag; cumulative exposure
  rejects(list(spec_model(M ~ lag2_A, var_type = "normal", mod_type = "mediator"), out),
          "its formula reads {lag2_A}")
  rejects(list(spec_model(M ~ cumA, var_type = "normal", mod_type = "mediator"), out),
          "its formula reads {cumA}")
  # the mediator's own recode reading the exposure (it runs on the
  # intervention's own data, before the swap)
  rejects(list(spec_model(M ~ AxV + V, var_type = "normal", mod_type = "mediator",
                          recode = recodes(AxV = A * V)), out),
          "its formula reads {AxV}")
  # an exposure lag carried inside the mediator's own recode
  rejects(list(spec_model(M ~ prevA, var_type = "normal", mod_type = "mediator",
                          recode = recodes(prevA = curA, curA = A)), out),
          "its formula reads {prevA}", init = recodes(curA = 0, prevA = 0), inr = NULL)
  # a self-referencing own recode
  rejects(list(spec_model(M ~ cum2, var_type = "normal", mod_type = "mediator",
                          recode = recodes(cum2 = cum2 + A)), out),
          "its formula reads {cum2}", init = recodes(cum2 = 0), inr = NULL)
  # a second in_recode entry transforming the lag (skipped by position, not name)
  rejects(list(spec_model(M ~ lag1_A, var_type = "normal", mod_type = "mediator"), out),
          "its formula reads {lag1_A}",
          init = recodes(lag1_A = 0), inr = recodes(lag1_A = A, lag1_A = 2 * lag1_A))
  # an out_recode copy: out_recode skips the first step, so not a lag
  rejects(list(spec_model(M ~ lagA, var_type = "normal", mod_type = "mediator"), out),
          "its formula reads {lagA}",
          init = recodes(lagA = 0), inr = NULL, outr = recodes(lagA = A))
  # a column created by another model's recode
  rejects(list(spec_model(L ~ A, var_type = "normal", mod_type = "covariate",
                          recode = recodes(A2 = A)),
               spec_model(M ~ A2, var_type = "normal", mod_type = "mediator"), out),
          "its formula reads {A2}", init = NULL, inr = NULL)
  # a subset reading the exposure or its lag (evaluated on the own-regime data)
  rejects(list(spec_model(M ~ V, var_type = "normal", mod_type = "mediator",
                          subset = A == 1), out),
          "its subset reads {A}")
  rejects(list(spec_model(M ~ V, var_type = "normal", mod_type = "mediator",
                          subset = lag1_A == 1), out),
          "its subset reads {lag1_A}")
})

testthat::test_that("check_natural_exposure_history rejects a mediator custom_sim it cannot inspect", {
  f   <- causalMed:::check_natural_exposure_history
  out <- spec_model(Y ~ A + M, var_type = "binary", mod_type = "outcome")
  csim <- function(fit, newdt) newdt[["cumA"]]

  # sim_value() hands custom_sim the whole swapped data set, so a formula scan
  # proves nothing about what it reads. With an exposure-derived column present
  # the pair is refused, even though the formula names only V.
  testthat::expect_error(
    f(list(spec_model(M ~ V, var_type = "custom", mod_type = "mediator",
                      custom_sim = csim), out),
      "A", recodes(cumA = 0), recodes(cumA = cumA + A), NULL),
    "supplies a custom_sim", fixed = TRUE)

  # Nothing is derived from the exposure here (the swap sets A and lag1_A), so
  # there is no stale column for it to read and the custom_sim is accepted.
  testthat::expect_no_error(
    f(list(spec_model(M ~ V + A + lag1_A, var_type = "custom",
                      mod_type = "mediator", custom_sim = csim), out),
      "A", recodes(lag1_A = 0), recodes(lag1_A = A), NULL))
})

# ---- Engine -----------------------------------------------------------------

# M depends on the exposure ONLY through its lag, so evaluating the mediator on
# the other regime's lag (Eq. 5) and keeping this regime's lag give far-apart
# answers with closed-form targets.
make_lag_med_data <- function(n = 3000, seed = 2026) {
  set.seed(seed)
  tt <- 0:2
  d <- data.table::CJ(id = seq_len(n), time = tt)
  data.table::setorderv(d, c("id", "time"))
  d[, V := rep(stats::rnorm(n), each = length(tt))]
  d[, A := stats::rbinom(.N, 1, 0.5)]
  d[, lag1_A := data.table::shift(A, fill = 0), by = id]
  d[, M := 2 * lag1_A + 0.5 * V + stats::rnorm(.N)]
  d[, Y := stats::rbinom(.N, 1, stats::plogis(-1 + 0.3 * A + 0.5 * M + 0.2 * V))]
  d[]
}

testthat::test_that("\"N\" evaluates the cross-world mediator on the other regime's exposure lag (Eq. 5)", {
  d <- make_lag_med_data()
  fit <- suppressWarnings(mediation(
    data = d, id_var = "id", base_vars = "V", exposure = "A", outcome = "Y",
    time_var = "time",
    models = list(
      spec_model(M ~ lag1_A + V, var_type = "normal", mod_type = "mediator"),
      spec_model(Y ~ A + M + V,  var_type = "binary", mod_type = "outcome")),
    init_recode = recodes(lag1_A = 0), in_recode = recodes(lag1_A = A),
    mediation_type = "N", mc_sample = 20000L, R = 1L, quiet = TRUE, seed = 7L))
  phi <- stats::setNames(fit$effect_size$Est, fit$effect_size$Intervention)

  # Targets at the last step (t = 2), A = 1:
  #   Eq. 5  : M = 0 + 0.5 V + e   (lag set to a' = 0)
  #   old    : M = 2 + 0.5 V + e   (lag kept at a = 1) -- what the current-
  #            exposure-only swap computed, NIE ~ 0
  set.seed(1)
  v <- stats::rnorm(2e5)
  e <- stats::rnorm(2e5)
  eq5 <- mean(stats::plogis(-0.7 + 0.5 * (0.5 * v + e) + 0.2 * v))
  old <- mean(stats::plogis(-0.7 + 0.5 * (2 + 0.5 * v + e) + 0.2 * v))
  testthat::expect_lt(abs(phi[["Phi10"]] - eq5), 0.04)
  testthat::expect_gt(abs(phi[["Phi10"]] - old), 0.15)
  nie <- fit$estimate$RD[fit$estimate$Effect == "Indirect effect"]
  testthat::expect_gt(nie, 0.15)   # truth ~0.22
})

testthat::test_that("the \"N\" swap leaves every input the exposure does not affect untouched", {
  # A mediator that does not depend on the exposure must be drawn identically
  # with and without the cross-world override. The first version of the
  # exposure-history fix re-ran the mediator's own recode on the swap copy,
  # which shifts an order-dependent chain like this one by a step.
  set.seed(5)
  n <- 400
  d <- data.table::data.table(id = seq_len(n), V = stats::rnorm(n),
                              A = stats::rbinom(n, 1, 0.5), L = stats::rnorm(n),
                              lag1_L = stats::rnorm(n), lag2_L = 0)
  d[, M := 0.5 * lag1_L + V + stats::rnorm(n)]
  fit <- causalMed:::fit_spec_models(
    list(spec_model(M ~ lag2_L + V, var_type = "normal", mod_type = "mediator",
                    recode = recodes(lag2_L = lag1_L, lag1_L = L))), d)
  draw <- function(interv) {
    set.seed(9)
    causalMed:::simulate_data(data.table::copy(d), exposure = "A", models = fit,
                              intervention = interv, mediation_type = "N")$M
  }
  own   <- draw(causalMed:::intervention_spec(1))                  # natural mediator
  cross <- draw(causalMed:::intervention_spec(1, list(M = 0)))     # cross-world swap
  testthat::expect_identical(cross, own)
})

testthat::test_that("mediation(\"N\") rejects exposure-derived inputs it cannot swap; \"I\" runs", {
  d <- make_lag_med_data(n = 300)
  d[, cumA := cumsum(A), by = id]
  d[, A_copy := A]
  base <- list(
    data = d, id_var = "id", base_vars = "V", exposure = "A", outcome = "Y",
    time_var = "time", mc_sample = 300L, R = 1L, quiet = TRUE, seed = 1L)
  m_out <- spec_model(Y ~ A + M + V, var_type = "binary", mod_type = "outcome")

  cum_args <- c(base, list(
    models = list(spec_model(M ~ cumA + V, var_type = "normal", mod_type = "mediator"), m_out),
    init_recode = recodes(cumA = 0), in_recode = recodes(cumA = cumA + A)))
  testthat::expect_error(do.call(mediation, c(cum_args, list(mediation_type = "N"))),
                         "its formula reads {cumA}", fixed = TRUE)
  testthat::expect_no_error(suppressWarnings(
    do.call(mediation, c(cum_args, list(mediation_type = "I")))))

  # A recode copying the current exposure: evaluated before the swap, so it
  # would keep the intervention's own exposure. Rejected; the formula term runs.
  testthat::expect_error(
    do.call(mediation, c(base, list(
      models = list(spec_model(M ~ A_copy + V, var_type = "normal", mod_type = "mediator",
                               recode = recodes(A_copy = A)), m_out),
      mediation_type = "N"))),
    "its formula reads {A_copy}", fixed = TRUE)
  testthat::expect_no_error(suppressWarnings(
    do.call(mediation, c(base, list(
      models = list(spec_model(M ~ A + V, var_type = "normal", mod_type = "mediator"), m_out),
      mediation_type = "N")))))
})

# ---- TMLE: model recodes ------------------------------------------------------

testthat::test_that(".tmle_check_model_recodes rejects model recodes derived from the exposure", {
  f <- causalMed:::.tmle_check_model_recodes
  testthat::expect_error(
    f(list(spec_model(M ~ AxV + V, var_type = "normal", mod_type = "mediator",
                      recode = recodes(AxV = A * V))), "A", NULL),
    "{AxV}", fixed = TRUE)
  # through a first-order lag, and transitively through another model's recode
  testthat::expect_error(
    f(list(spec_model(L ~ V, var_type = "normal", mod_type = "covariate",
                      recode = recodes(LA = lag_A * 2)),
           spec_model(Y ~ LA2, var_type = "binary", mod_type = "outcome",
                      recode = recodes(LA2 = LA + 1))),
      "A", recodes(lag_A = A)),
    "{LA, LA2}", fixed = TRUE)
  # recodes that do not touch the exposure pass
  testthat::expect_no_error(
    f(list(spec_model(M ~ V_sq, var_type = "normal", mod_type = "mediator",
                      recode = recodes(V_sq = V^2))), "A", recodes(lag_A = A)))
})

testthat::test_that("mediation(estimator = \"tmle\") refuses a model recode that reads the exposure", {
  d <- make_lag_med_data(n = 300)
  d[, AxV := A * V]
  testthat::expect_error(suppressWarnings(mediation(
    data = d, id_var = "id", base_vars = "V", exposure = "A", outcome = "Y",
    time_var = "time",
    models = list(
      spec_model(A ~ V, var_type = "binary", mod_type = "exposure"),
      spec_model(M ~ AxV + V, var_type = "normal", mod_type = "mediator",
                 recode = recodes(AxV = A * V)),
      spec_model(Y ~ A + M + V, var_type = "binary", mod_type = "outcome")),
    mediation_type = "N", estimator = "tmle", quiet = TRUE)),
    "{AxV}", fixed = TRUE)
})
