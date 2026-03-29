# Unit tests for mediation()
# Only point estimates (R = 1); bootstrap CIs should be computed locally, not in CI.
# Requires: testthat (>= 3.0.0)


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
  data(nonsurvivaldata)

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
        R              = 0L,
        quiet          = TRUE,
        seed           = 42L
      )
    )

    est      <- fit$estimate
    total    <- est$RD[est$Effect == "Total effect"]
    direct   <- est$RD[est$Effect == "Direct effect"]
    indirect <- est$RD[est$Effect == "Indirect effect"]

    # Exact algebraic identity: must hold to floating-point precision
    testthat::expect_equal(
      direct + indirect, total,
      tolerance = 1e-10,
      label = paste0("Direct + Indirect = Total [type = ", mtype, "]")
    )

    # Mediation proportion = Indirect / Total * 100
    med_prop      <- est$RD[est$Effect == "Mediation Proportion"]
    expected_prop <- if (isTRUE(abs(total) < 1e-10)) NA_real_ else indirect / total * 100
    testthat::expect_equal(
      med_prop, expected_prop,
      tolerance = 1e-10,
      label = paste0("Mediation Proportion [type = ", mtype, "]")
    )
  }
})

