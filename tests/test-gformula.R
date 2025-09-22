# Unit tests for gformula()
# Only point estimates (R = 1); bootstrap CIs should be computed locally, not in CI.
# Requires: testthat (>= 3.0.0)
testthat::local_edition(3)
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

testthat::test_that("Fixed intervention", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("causalMed")
  
  suppressPackageStartupMessages(library(causalMed))
  suppressPackageStartupMessages(library(gfoRmula))
  
  data("nonsurvivaldata", package = "causalMed")
  dat <- data.table::as.data.table(nonsurvivaldata)
  dat[, time := as.integer(time)]
  data.table::setorder(dat, id, time)

  time_points <- max(dat$time, na.rm = TRUE) + 1L
  
  # ----- recodes for lags (causalMed side) ---------------------------------
  init_recodes <- c("lag1_A = 0", "lag1_L1 = 0", "lag1_L2 = 0")
  in_recodes   <- c("lag1_A = A", "lag1_L1 = L1", "lag1_L2 = L2")
  
  # ----- models (temporal order) -------------------------------------------
  m_L1 <- spec_model(
    L1 ~ lag1_A + lag1_L1 + V + time,
    var_type = "normal",   mod_type = "covariate"
  )
  m_L2 <- spec_model(
    L2 ~ lag1_A + lag1_L2 + V + time,
    var_type = "binomial", mod_type = "covariate"
  )
  m_A <- spec_model(
    A  ~ lag1_A + L1 + L2 + V + time,
    var_type = "binomial", mod_type = "exposure"
  )
  m_Y <- spec_model(
    Y_bin  ~ A + L1 + L2 ,
    var_type = "binomial", mod_type = "outcome"
  )
  models <- list(m_L1, m_L2, m_A, m_Y)
  
  # ----- causalMed::gformula (natural vs always=1) --------------------------
  mc_nsim <- 5000L  # Monte Carlo size; balance accuracy vs test time
  fit_cm <- testthat::expect_no_error(
    causalMed::gformula(
      data          = dat,
      id_var        = "id",
      time_var      = "time",
      base_vars     = "V",
      exposure      = "A",
      models        = models,
      intervention  = list(natural = NULL, always = 1),
      init_recode   = init_recodes,
      in_recode     = in_recodes,
      out_recode    = NULL,
      mc_sample     = mc_nsim,
      R             = 1,           # point estimates only
      return_data   = FALSE,
      return_fitted = FALSE,
      quiet         = TRUE,
      seed          = 20250915
    )
  )
  es_cm <- data.table::as.data.table(fit_cm$effect_size)
  means_cm <- as.numeric(es_cm$Est[1:2])
  
  # ----- gfoRmula::gformula (natural vs A=1 static) -------------------------
  dat2<-dat
  dat2$time<-as.integer(dat2$time)
  res_gf <- testthat::expect_no_error(
    gfoRmula::gformula(
      obs_data      = dat2,
      id            = "id",
      time_name     = "time",
      time_points   = 5,
      covnames      = c('L1', 'L2', 'A'),
      covtypes      = c("normal", "binary", "binary"),
      basecovs      = c("V"),
      covparams     = list(covmodels = c(
        L1 ~ lag1_A + lag1_L1 + V + time,
        L2 ~ lag1_A + lag1_L2 + V + time,
        A  ~ lag1_A + L1 + L2 + V + time
      )),
      ymodel        = Y_bin ~ A + L1 + L2 ,
      outcome_name  = "Y_bin",
      outcome_type  = "binary_eof",
      histories     = c(gfoRmula::lagged),
      histvars      = list(c('A', 'L1', 'L2')),
      intervention1.A = list(gfoRmula::static, rep(1, 5)),  
      int_descript  = "Always treat",
      nsimul        = 50000,
      seed          = 20250915
    )
  )
  
  gfm <- res_gf$result[["g-form mean"]]
  means_gf <- as.numeric(gfm[1:2])
  
  # ----- Compare with tolerance --------------------------------------------
  tol <- 0.01
  testthat::expect_equal(means_cm, means_gf, tolerance = tol)
})


testthat::test_that("Predefined intervention", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("causalMed")
  
  suppressPackageStartupMessages(library(causalMed))
  suppressPackageStartupMessages(library(gfoRmula))
  
  data("nonsurvivaldata", package = "causalMed")
  dat <- data.table::as.data.table(nonsurvivaldata)
  dat[, time := as.integer(time)]
  data.table::setorder(dat, id, time)
  
  time_points <- max(dat$time, na.rm = TRUE) + 1L
  
  # ----- recodes for lags (causalMed side) ---------------------------------
  init_recodes <- c("lag1_A = 0", "lag1_L1 = 0", "lag1_L2 = 0")
  in_recodes   <- c("lag1_A = A", "lag1_L1 = L1", "lag1_L2 = L2")
  
  # ----- models (temporal order) -------------------------------------------
  m_L1 <- spec_model(
    L1 ~ lag1_A + lag1_L1 + V + time,
    var_type = "normal",   mod_type = "covariate"
  )
  m_L2 <- spec_model(
    L2 ~ lag1_A + lag1_L2 + V + time,
    var_type = "binomial", mod_type = "covariate"
  )
  m_A <- spec_model(
    A  ~ lag1_A + L1 + L2 + V + time,
    var_type = "binomial", mod_type = "exposure"
  )
  m_Y <- spec_model(
    Y_bin  ~ A + L1 + L2 ,
    var_type = "binomial", mod_type = "outcome"
  )
  models <- list(m_L1, m_L2, m_A, m_Y)
  
  # ----- causalMed::gformula (natural vs always=1) --------------------------
  mc_nsim <- 5000L  # Monte Carlo size; balance accuracy vs test time
  fit_cm <- testthat::expect_no_error(
    causalMed::gformula(
      data          = dat,
      id_var        = "id",
      time_var      = "time",
      base_vars     = "V",
      exposure      = "A",
      models        = models,
      intervention  = list(Predefined=c(0,0,1,1,1)),
      init_recode   = init_recodes,
      in_recode     = in_recodes,
      out_recode    = NULL,
      mc_sample     = mc_nsim,
      R             = 1,           # point estimates only
      return_data   = FALSE,
      return_fitted = FALSE,
      quiet         = TRUE,
      seed          = 20250915
    )
  )
  es_cm <- data.table::as.data.table(fit_cm$effect_size)
  means_cm <- as.numeric(es_cm$Est[1:2])
  
  # ----- gfoRmula::gformula (natural vs A=1 static) -------------------------
  dat2<-dat
  dat2$time<-as.integer(dat2$time)
  res_gf <- testthat::expect_no_error(
    gfoRmula::gformula(
      obs_data      = dat2,
      id            = "id",
      time_name     = "time",
      time_points   = 5,
      covnames      = c('L1', 'L2', 'A'),
      covtypes      = c("normal", "binary", "binary"),
      basecovs      = c("V"),
      covparams     = list(covmodels = c(
        L1 ~ lag1_A + lag1_L1 + V + time,
        L2 ~ lag1_A + lag1_L2 + V + time,
        A  ~ lag1_A + L1 + L2 + V + time
      )),
      ymodel        = Y_bin ~ A + L1 + L2 ,
      outcome_name  = "Y_bin",
      outcome_type  = "binary_eof",
      histories     = c(gfoRmula::lagged),
      histvars      = list(c('A', 'L1', 'L2')),
      intervention1.A = list(gfoRmula::static, c(0,0,1,1,1)),  
      int_descript  = "Always treat",
      nsimul        = 50000,
      seed          = 20250915
    )
  )
  
  gfm <- res_gf$result[["g-form mean"]]
  means_gf <- as.numeric(gfm[1:2])
  
  # ----- Compare with tolerance --------------------------------------------
  tol <- 0.01
  testthat::expect_equal(means_cm, means_gf, tolerance = tol)
})



testthat::test_that("Dynamic intervention: causalMed vs gfoRmula", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("causalMed")
  testthat::skip_if_not_installed("gfoRmula")
  
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(causalMed))
  suppressPackageStartupMessages(library(gfoRmula))
  
  data("binary_eofdata", package = "gfoRmula")
  dat <- as.data.table(binary_eofdata)
  setorder(dat, id_num, time)
  time_points <- max(dat$time, na.rm = TRUE) + 1L
  
  binary_eofdata[, `:=`(
    lag1_treat = shift(treat, 1L, type = "lag", fill=0), 
    lag1_cov1 = shift(cov1, 1L, type = "lag", fill=0),
    lag1_cov2 = shift(cov2, 1L, type = "lag", fill=0)
  ), by = id_num]
  
  # ----- recodes for lags (causalMed side) ---------------------------------
  init_recodes <- c("lag1_treat = 0", "lag1_cov1 = 0", "lag1_cov2 = 0")
  in_recodes   <- c("lag1_treat = treat", "lag1_cov1 = cov1", "lag1_cov2 = cov2")
  
  # ----- models (temporal order) -------------------------------------------
  m_cov1 <- spec_model(
    cov1 ~ lag1_treat + lag1_cov1 + cov3 + time,
    var_type = "binomial", mod_type = "covariate"
  )
  m_cov2 <- spec_model(
    cov2 ~ lag1_treat + lag1_cov2 + cov3 + time,
    var_type = "normal",   mod_type = "covariate"
  )
  m_A <- spec_model(
    treat ~ lag1_treat + cov1 + cov2 + cov3 + time,
    var_type = "normal", mod_type = "exposure"
  )
  m_Y <- spec_model(
    outcome ~ treat + cov1 + cov2,
    var_type = "binomial", mod_type = "outcome"
  )
  models <- list(m_cov1, m_cov2, m_A, m_Y)
  
  set.seed(20250915)
  mc_nsim_cm <- 5000L   
  mc_nsim_gf <- 5000L  
  
  ## ---- causalMed::gformula----
  fit_cm <- testthat::expect_no_error(
    causalMed::gformula(
      data          = dat,
      id_var        = "id_num",
      time_var      = "time",
      base_vars     = "cov3",
      exposure      = "treat",
      models        = models,
      intervention  = list(dynamic = "treat>0"),
      init_recode   = init_recodes,
      in_recode     = in_recodes,
      out_recode    = NULL,
      mc_sample     = mc_nsim_cm,
      R             = 1,          # point estimates only
      return_data   = F,
      return_fitted = F,
      quiet         = TRUE,
      seed          = 20250915
    )
  )
  es_cm <- as.data.table(fit_cm$effect_size)
  means_cm <- as.numeric(es_cm$Est[1:2])
  
  # ----- gfoRmula::gformula (natural vs A=1 static) -------------------------
  treat_gt0 <- function(newdf, pool, intvar, intvals, time_name, t){
    newdf[, (intvar) := as.numeric(newdf[[intvar]] > 0)]
  }
  res_gf <- testthat::expect_no_error(
    gfoRmula::gformula(
      obs_data      = dat,
      id            = "id_num",
      time_name     = "time",
      time_points   = time_points,
      covnames      = c("cov1", "cov2", "treat"),
      covtypes      = c("binary", "normal", "normal"),
      basecovs      = "cov3",
      covparams     = list(covmodels = c(
        cov1  ~ lag1_treat + lag1_cov1 + cov3 + time,
        cov2  ~ lag1_treat + lag1_cov2 + cov3 + time,
        treat ~ lag1_treat + cov1 + cov2 + cov3 + time
      )),
      ymodel        = outcome ~ treat + cov1 + cov2,
      outcome_name  = "outcome",
      outcome_type  = "binary_eof",
      histories     = c(gfoRmula::lagged),
      histvars      = list(c("treat","cov1","cov2")),
      intervention1.treat = list(treat_gt0),   
      int_descript = "Treat=1 if natural treat>0; else 0",
      nsimul        = mc_nsim_gf,
      seed          = 20250915
    )
  )
  res_tab <- as.data.table(res_gf$result)
  gform_col <- grep("^g-form", names(res_tab), value = TRUE)[1]
  means_gf <- as.numeric(res_tab[[gform_col]][1:2])
  
  # ----- Compare with tolerance --------------------------------------------
  tol <- 0.01
  testthat::expect_equal(means_cm, means_gf, tolerance = tol)
})
