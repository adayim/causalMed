
testthat::test_that("mediation runs without error on nonsurvivaldata", {
  testthat::skip_on_cran()

  data("nonsurvivaldata", package = "causalMed")

  models<-list(
    modA<-spec_model(A  ~ V + lag1_L1 + lag1_L2 + lag1_A + time,var_type  = "binary",mod_type = "exposure"),
    modL1<-spec_model(L1 ~ V + A + lag1_L1 + time,var_type  = "normal",mod_type = "covariate"),
    modL2<-spec_model(L2 ~ V + A + lag1_L2 + time,var_type  = "binary",mod_type = "covariate"),
    modM<-spec_model(M  ~ V + A + L1 + L2 + lag1_M + time,var_type  = "normal",mod_type = "mediator"),
    modY<-spec_model(Y_bin  ~ V + A + M + L1 + L2 + A:M,var_type  = "binary",mod_type = "outcome")#Y_cont同
  )

  testthat::expect_no_error(
    mediation(
      data         = nonsurvivaldata,
      id_var       = "id",
      time_var     = "time",
      base_vars    = "V",
      exposure     = "A",
      outcome      = "Y_bin",
      models       = models,
      init_recode  = recodes(lag1_A=0,lag1_L1=0,lag1_M=0,lag1_L2=0),
      in_recode    = recodes(lag1_A=A,lag1_L1=L1,lag1_M=M,lag1_L2=L2),
      out_recode   = NULL,
      mediation_type = "I",
      mc_sample    = 90000,
      R            = 1,
      return_data  = FALSE,
      return_fitted = FALSE
    )
  )
})



testthat::test_that("Mediation decomposition: Direct + Indirect = Total (algebraic identity)", {
  data("nonsurvivaldata", package = "causalMed") 

  # Minimal models: A (exposure), L2 (confounder), L1 (mediator), Y_bin (outcome)
  # Ordering: A -> L2 -> M(L1) -> Y
  m_L2  <- spec_model(L2    ~ A + V + time,         var_type = "binary", mod_type = "covariate")
  m_med <- spec_model(L1    ~ A + V + L2 + time,    var_type = "normal", mod_type = "mediator")
  m_Y   <- spec_model(Y_bin ~ A + L1 + L2 + V,      var_type = "binary", mod_type = "outcome")
  models <- list(m_L2, m_med, m_Y)

  for (mtype in c("N", "I")) {
    fit <- suppressWarnings(
      mediation(
        data           = nonsurvivaldata,
        id_var         = "id",
        base_vars      = "V",
        exposure       = "A",
        outcome        = "Y_bin",
        time_var       = "time",
        models         = models,
        mediation_type = mtype,
        mc_sample      = 500L,
        R              = 1L,
        quiet          = TRUE,
        seed           = 42L
      )
    )

    est      <- fit$estimate
    total    <- est$RD[est$Effect == "Total effect"]
    direct   <- est$RD[est$Effect == "Direct effect"]
    indirect <- est$RD[est$Effect == "Indirect effect"]
    resid    <- est$RD[est$Effect == "TE - (Direct + Indirect)"]

    if (mtype == "N") {
      # Natural effects (Zheng & van der Laan 2017): the decomposition sums
      # exactly to the total effect; no residual row.
      testthat::expect_false("TE - (Direct + Indirect)" %in% est$Effect)
      testthat::expect_equal(
        direct + indirect, total, tolerance = 1e-10,
        label = "Direct + Indirect = Total [type = N]"
      )
    } else {
      # Interventional effects (Lin 2017; VanderWeele & Tchetgen 2017; Yamamuro
      # 2021): Direct + Indirect = the interventional overall effect, and the
      # natural total effect is recovered only after adding the residual.
      testthat::expect_length(resid, 1L)
      testthat::expect_equal(
        direct + indirect + resid, total, tolerance = 1e-10,
        label = "Direct + Indirect + residual = Total [type = I]"
      )
    }

    # Additive mediation proportion = (Total - Direct) / Total * 100
    # (= (sum IIE + residual) / Total; reduces to Indirect/Total when type = N).
    med_prop      <- est$RD[est$Effect == "Mediation Proportion"]
    expected_prop <- if (isTRUE(abs(total) < 1e-10)) NA_real_ else (total - direct) / total * 100
    testthat::expect_equal(
      med_prop, expected_prop,
      tolerance = 1e-10,
      label = paste0("Mediation Proportion [type = ", mtype, "]")
    )
  }
})


testthat::test_that("Multi-mediator (N=2, type=I): sum(IIE_k) + IDE = TE (Yamamuro 2021)", {
  data("nonsurvivaldata", package = "causalMed")

  # Two mediators in temporal order: A -> L1 -> L2 -> Y_bin
  m_M1 <- spec_model(L1    ~ V + A + time,             var_type = "normal", mod_type = "mediator")
  m_M2 <- spec_model(L2    ~ V + A + L1 + time,        var_type = "binary", mod_type = "mediator")
  m_Y  <- spec_model(Y_bin ~ V + A + L1 + L2,          var_type = "binary", mod_type = "outcome")

  fit <- suppressWarnings(
    mediation(
      data           = nonsurvivaldata,
      id_var         = "id",
      base_vars      = "V",
      exposure       = "A",
      outcome        = "Y_bin",
      time_var       = "time",
      models         = list(m_M1, m_M2, m_Y),
      mediation_type = "I",
      n_vw           = 1L,           # speed: single permutation
      mc_sample      = 1000L,
      R              = 1L,
      quiet          = TRUE,
      seed           = 7L
    )
  )

  est <- fit$estimate

  # Sequential interventional decomposition (Yamamuro 2021): the per-mediator
  # IIEs and the IDE sum to the interventional OVERALL effect; the natural total
  # effect is recovered after adding the mediated-interaction residual.
  total    <- est$RD[est$Effect == "Total effect"]
  direct   <- est$RD[est$Effect == "Direct effect"]
  resid    <- est$RD[est$Effect == "TE - (Direct + Indirect)"]
  iie_rd   <- est$RD[grepl("^Indirect effect", est$Effect)]

  testthat::expect_length(iie_rd, 2L)   # one IIE row per mediator
  testthat::expect_length(resid, 1L)
  testthat::expect_equal(
    direct + sum(iie_rd) + resid, total,
    tolerance = 1e-10,
    label = "sum(IIE_k) + IDE + residual = TE [N=2, I]"
  )

  # Arms: natural-course nat0/nat1 (for TE) plus permuted-pool interventional
  # interventions Phi00, Phi10, Phi1_1, Phi11.
  testthat::expect_setequal(
    fit$effect_size$Intervention,
    c("nat0", "nat1", "Phi00", "Phi10", "Phi1_1", "Phi11")
  )

  # Multiplicative PM row is present
  testthat::expect_true("Mediation Proportion (multiplicative)" %in% est$Effect)
})


testthat::test_that("Multi-mediator with mediation_type='N' is rejected", {
  data("nonsurvivaldata", package = "causalMed")
  m_M1 <- spec_model(L1    ~ V + A + time,      var_type = "normal", mod_type = "mediator")
  m_M2 <- spec_model(L2    ~ V + A + L1 + time, var_type = "binary", mod_type = "mediator")
  m_Y  <- spec_model(Y_bin ~ V + A + L1 + L2,   var_type = "binary", mod_type = "outcome")

  testthat::expect_error(
    suppressWarnings(mediation(
      data           = nonsurvivaldata,
      id_var         = "id", base_vars = "V",
      exposure       = "A",  outcome   = "Y_bin",
      time_var       = "time",
      models         = list(m_M1, m_M2, m_Y),
      mediation_type = "N",
      mc_sample      = 100L,
      R              = 1L,
      quiet          = TRUE
    )),
    regexp = "Multiple mediators"
  )
})



testthat::test_that("mediation retains per-replicate bootstrap estimates and prints new fields", {
  data("nonsurvivaldata", package = "causalMed")

  m_L2  <- spec_model(L2    ~ A + V + time,      var_type = "binary", mod_type = "covariate")
  m_med <- spec_model(L1    ~ A + V + L2 + time, var_type = "normal", mod_type = "mediator")
  m_Y   <- spec_model(Y_bin ~ A + L1 + L2 + V,   var_type = "binary", mod_type = "outcome")
  models <- list(m_L2, m_med, m_Y)

  fit <- suppressWarnings(
    mediation(
      data           = nonsurvivaldata,
      id_var         = "id",
      base_vars      = "V",
      exposure       = "A",
      outcome        = "Y_bin",
      time_var       = "time",
      models         = models,
      mediation_type = "I",
      mc_sample      = 500L,
      R              = 4L,
      quiet          = TRUE,
      seed           = 42L
    )
  )

  # Interventions: 4 replicates x 5 (nat0, nat1, Phi00, Phi10, Phi11)
  bi <- fit$boot_estimates$interventions
  testthat::expect_true(all(c("replicate", "Intervention", "Est") %in% names(bi)))
  testthat::expect_equal(nrow(bi), 4L * 5L)

  # Effects: 4 replicates x (Indirect, Direct, Total, residual + 2 PM rows) = 4 x 6
  be <- fit$boot_estimates$effects
  testthat::expect_true(all(c("replicate", "Effect", "RD", "RR") %in% names(be)))
  testthat::expect_equal(nrow(be), 4L * 6L)
  testthat::expect_true(all(c("Mediation Proportion",
                              "Mediation Proportion (multiplicative)") %in% be$Effect))

  # Data summary and observed benchmark
  testthat::expect_equal(fit$data_summary$n_id, length(unique(nonsurvivaldata$id)))
  testthat::expect_true(is.finite(fit$observed$value))

  # print: mediator line, n_vw, observed benchmark, CI footnote,
  # correctly-scaled decomposition CI labels
  out <- paste(capture.output(print(fit)), collapse = "\n")
  testthat::expect_match(out, "Mediator(s)  : L1", fixed = TRUE)
  testthat::expect_match(out, "Outcome      : Y_bin")
  testthat::expect_match(out, "n_vw")
  testthat::expect_match(out, "Observed (nonparametric)", fixed = TRUE)
  testthat::expect_match(out, "RD 2.5%", fixed = TRUE)
  testthat::expect_match(out, "RR 2.5%", fixed = TRUE)
  testthat::expect_match(out, "95% CIs", fixed = TRUE)
})



testthat::test_that("natural-effect intermediate-confounder caveat is retained and printed", {
  data("nonsurvivaldata", package = "causalMed")

  # L2's covariate model includes the exposure A on the RHS -> exposure-affected
  # (intermediate) confounder, so the natural effects are not point-identified.
  m_L2  <- spec_model(L2    ~ A + V + time,      var_type = "binary", mod_type = "covariate")
  m_med <- spec_model(L1    ~ A + V + L2 + time, var_type = "normal", mod_type = "mediator")
  m_Y   <- spec_model(Y_bin ~ A + L1 + L2 + V,   var_type = "binary", mod_type = "outcome")
  models <- list(m_L2, m_med, m_Y)

  fit <- suppressWarnings(
    mediation(
      data           = nonsurvivaldata,
      id_var         = "id",
      base_vars      = "V",
      exposure       = "A",
      outcome        = "Y_bin",
      time_var       = "time",
      models         = models,
      mediation_type = "N",
      mc_sample      = 500L,
      R              = 1L,
      quiet          = TRUE,
      seed           = 42L
    )
  )

  # The offending confounder is retained on the object ...
  testthat::expect_true("L2" %in% fit$intermediate_confounders)
  # ... and re-surfaced by print() in a short form.
  out <- paste(capture.output(print(fit)), collapse = "\n")
  testthat::expect_match(out, "Identifiability", fixed = TRUE)
  testthat::expect_match(out, "not point-identified", ignore.case = TRUE)

  # Interventional analyses do not carry the caveat.
  fit_i <- suppressWarnings(
    mediation(
      data = nonsurvivaldata, id_var = "id", base_vars = "V", exposure = "A",
      outcome = "Y_bin", time_var = "time", models = models,
      mediation_type = "I", mc_sample = 500L, R = 1L, quiet = TRUE, seed = 42L
    )
  )
  out_i <- paste(capture.output(print(fit_i)), collapse = "\n")
  testthat::expect_false(grepl("Identifiability", out_i, fixed = TRUE))
})
