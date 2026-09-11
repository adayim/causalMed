
test_that("gformula throws errors for invalid data input", {
  # Mock valid inputs to use as a base
  data(nonsurvivaldata)
  
  # Models need to be spec_model
  m1 <- spec_model(A ~ V, var_type = "binary", mod_type = "exposure")
  m2 <- spec_model(Y_bin ~ A + L1 + V, var_type = "binary", mod_type = "outcome")
  models <- list(m1, m2)
  
  # data not dataframe
  expect_error(
    gformula(data = "not a dataframe", 
             id_var = "id", base_vars = "V", exposure = "A", time_var = "time", models = models),
    "Data does not exist or not a data.frame"
  )
  
  # time_var not numeric - nonsurvivaldata time is int (numeric), let's make a broken one
  bad_data <- copy(nonsurvivaldata)
  bad_data$time_char <- as.character(bad_data$time)
  expect_error(
    gformula(data = bad_data, 
             id_var = "id", base_vars = "V", exposure = "A", time_var = "time_char", models = models),
    "The time variable must be numeric"
  )
  
  # check_var_in failures
  expect_error(
    gformula(data = nonsurvivaldata, 
             id_var = "bad_id", base_vars = "V", exposure = "A", time_var = "time", models = models),
    "The following variables cannot be found in the data: bad_id"
  )
})

test_that("gformula throws errors for invalid model specifications", {
  data(nonsurvivaldata)
  m1 <- spec_model(A ~ V, var_type = "binary", mod_type = "exposure")
  m_out <- spec_model(Y_bin ~ A + L1 + V, var_type = "binary", mod_type = "outcome")
  
  # Models not list
  expect_error(
    gformula(data = nonsurvivaldata, 
             id_var = "id", base_vars = "V", exposure = "A", time_var = "time", models = "not list"),
    "Models must be provided as a list"
  )
  
  # Multiple outcome models
  m_out2 <- spec_model(Y_bin ~ A, var_type = "binary", mod_type = "outcome")
  expect_error(
     gformula(data = nonsurvivaldata, 
             id_var = "id", base_vars = "V", exposure = "A", time_var = "time", models = list(m1, m_out, m_out2)),
     "Only one outcome model or survival model allowed"
  )
  
  # No outcome model
  expect_error(
     gformula(data = nonsurvivaldata, 
             id_var = "id", base_vars = "V", exposure = "A", time_var = "time", models = list(m1)),
     "Outcome model or survival model must be defined"
  )
  
  # Multiple exposure models
  m2 <- spec_model(A ~ L1, var_type = "binary", mod_type = "exposure")
  expect_error(
    gformula(data = nonsurvivaldata, 
             id_var = "id", base_vars = "V", exposure = "A", time_var = "time", models = list(m1, m2, m_out)),
    "Only one exposure model is allowed"
  )
  
  # Exposure variable mismatch
  m3 <- spec_model(L1 ~ V, var_type = "binary", mod_type = "exposure") # incorrectly labeled check
  # Wait, spec_model doesn't enforce variable name matching, but gformula check_error does:
  # "The given exposure variable was different between exposure model in `models`."
  expect_error(
    gformula(data = nonsurvivaldata, 
             id_var = "id", base_vars = "V", exposure = "A", time_var = "time", models = list(m3, m_out)),
    "The given exposure variable was different between exposure model"
  )
  
})

test_that("gformula throws errors for invalid intervention inputs", {
  data(nonsurvivaldata)
  m1 <- spec_model(A ~ V, var_type = "binary", mod_type = "exposure")
  m2 <- spec_model(L1 ~ A + V, var_type="normal", mod_type="covariate")
  m3 <- spec_model(Y_bin ~ A + L1 + V, var_type = "binary", mod_type = "outcome")
  models <- list(m1, m2, m3)

  # Check with mediator in models but intervention not NULL
  m_med <- spec_model(M ~ A + L1, var_type="normal", mod_type="mediator")
  models_med <- list(m1, m2, m_med, m3)
  
  expect_error(
     gformula(data = nonsurvivaldata, 
             id_var = "id", base_vars = "V", exposure = "A", time_var = "time", models = models_med,
             intervention = list(always = 1)),
     "You cannot specify intervention with mediator at the same time"
  )
  
  # Intervention not a list
  expect_error(
    gformula(data = nonsurvivaldata, models = models, id_var="id", base_vars="V", exposure="A", time_var="time",
             intervention = "not a list"),
    "Intervention must be a list object"
  )
  
  # Unknown ref_int
  expect_error(
    gformula(data = nonsurvivaldata, models = models, id_var="id", base_vars="V", exposure="A", time_var="time",
             intervention = list(always=1), ref_int = "never"),
    "`ref_int` must be included in the names of intervention"
  )
  
})

test_that("gformula throws error for invalid recode parameters", {
  data(nonsurvivaldata)
  m1 <- spec_model(A ~ V, var_type = "binary", mod_type = "exposure")
  m3 <- spec_model(Y_bin ~ A + L1 + V, var_type = "binary", mod_type = "outcome")
  models <- list(m1, m3)
  
  # Passed as character
  expect_error(
    gformula(data = nonsurvivaldata, 
             id_var = "id", base_vars = "V", exposure = "A", time_var = "time", models = models,
             init_recode = "A = 1"),
    "Invalid input for 'init_recode'"
  )
  
  # Passed as list (but not causalMed_recodes class)
  expect_error(
    gformula(data = nonsurvivaldata,
             id_var = "id", base_vars = "V", exposure = "A", time_var = "time", models = models,
             in_recode = list(A = 1)),
    "Invalid input for 'in_recode'"
  )

})

test_that("mediation validates exposure_regime / reference_regime", {
  data(nonsurvivaldata)   # 5 distinct time points
  m_L2  <- spec_model(L2    ~ A + V + time,      var_type = "binary", mod_type = "covariate")
  m_med <- spec_model(L1    ~ A + V + L2 + time, var_type = "normal", mod_type = "mediator")
  m_Y   <- spec_model(Y_bin ~ A + L1 + L2 + V,   var_type = "binary", mod_type = "outcome")
  models <- list(m_L2, m_med, m_Y)
  med <- function(...) mediation(
    data = nonsurvivaldata, id_var = "id", base_vars = "V", exposure = "A",
    outcome = "Y_bin", time_var = "time", models = models,
    mc_sample = 100L, R = 1L, quiet = TRUE, ...)

  expect_error(med(exposure_regime = c(0, 1, 1)),
               "`exposure_regime` must have length 1 or 5 (one value per distinct time point); got length 3",
               fixed = TRUE)
  expect_error(med(reference_regime = c(0, 0, 2, 0, 0)),
               "`reference_regime` must take values in {0, 1}; got {2}", fixed = TRUE)
  expect_error(med(exposure_regime = "1"),
               "`exposure_regime` must be a numeric or logical vector", fixed = TRUE)
  expect_error(med(exposure_regime = c(0, 0, 0, 0, 0)),
               "`exposure_regime` and `reference_regime` are identical", fixed = TRUE)
  expect_error(med(exposure_regime = 1, reference_regime = TRUE),
               "are identical", fixed = TRUE)
  # TMLE: defaults only. suppressWarnings(): the "N" identifiability warning
  # (L2 ~ A) fires before the TMLE block is reached.
  m_A <- spec_model(A ~ V + time, var_type = "binary", mod_type = "exposure")
  expect_error(suppressWarnings(
    mediation(data = nonsurvivaldata, id_var = "id", base_vars = "V", exposure = "A",
              outcome = "Y_bin", time_var = "time", models = c(list(m_A), models),
              mediation_type = "N", estimator = "tmle", exposure_regime = c(0, 1, 1, 0, 0),
              quiet = TRUE)),
    "estimator = 'tmle' supports only the default regimes", fixed = TRUE)
})
