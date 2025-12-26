# Unit tests for mediation()
# Only point estimates (R = 1); bootstrap CIs should be computed locally, not in CI.
# Requires: testthat (>= 3.0.0)


testthat::test_that("mediation runs without error on nonsurvivaldata", {
  testthat::skip_on_cran()

  data("nonsurvivaldata", package = "causalMed")

  models<-list(
    modA<-spec_model(A  ~ V + lag1_L1 + lag1_L2 + lag1_A + time,var_type  = "binomial",mod_type = "exposure"),
    modL1<-spec_model(L1 ~ V + A + lag1_L1 + time,var_type  = "normal",mod_type = "covariate"),
    modL2<-spec_model(L2 ~ V + A + lag1_L2 + time,var_type  = "binomial",mod_type = "covariate"),
    modM<-spec_model(M  ~ V + A + L1 + L2 + lag1_M + time,var_type  = "normal",mod_type = "mediator"),
    modY<-spec_model(Y_bin  ~ V + A + M + L1 + L2 + A:M,var_type  = "binomial",mod_type = "outcome")#Y_cont同
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
