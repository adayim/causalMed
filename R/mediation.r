
#' G-formula Analysis
#'
#' @description Mediation analysis for time varying mediator, estimation based
#'  on g-formula. Output contains total effect, #' natural direct effect and natural
#'   indirect effect for mediation or regular g-formula. data.frame will be returned.
#'
#' @note Not that final outcome must be the same for in all rows per subject.
#'   If the dataset is  survival settings, records after the interested outcome
#'    must be deleted. The function it self do some data manipulation internally.
#'    Please prepare the data as longitudinal format.
#'
#' @inheritParams gformula
#'
#' @param mediation_type Type of the mediation effect, natural effect (\code{"N"}) or interventional effect (\code{"I"}).
#'
#' @references
#' Lin, S. H., Young, J. G., Logan, R., & VanderWeele, T. J. (2017). Mediation analysis for a survival outcome with time-varying exposures, mediators, and confounders. \emph{Statistics in medicine}, 36(26), 4153-4166. DOI:10.1002/sim.7426
#' Zheng, W., & van der Laan, M. (2017). Longitudinal mediation analysis with time-varying mediators and exposures, with application to survival outcomes. \emph{Journal of causal inference}, 5(2). DOI:10.1515/jci-2016-0006
#'
#'
#' TODO: weights, time varying intervention
#'
#' @export
#'
#'
mediation <- function(data,
                      id_var,
                      base_vars,
                      exposure,
                      outcome,
                      time_var,
                      models,
                      init_recode = NULL,
                      in_recode = NULL,
                      out_recode = NULL,
                      mc_sample = nrow(data),
                      mediation_type = c("N", "I"),
                      return_fitted = FALSE,
                      R = 500,
                      ncores = 1L) {
  tpcall <- match.call()
  mediation_type <- match.arg(mediation_type)

  # Calculate mediation effect
  risk_calc <- function(data_list) {
    phi_11 <- mean(data_list$always[["Pred_Y"]])
    phi_00 <- mean(data_list$never[["Pred_Y"]])
    phi_10 <- mean(data_list$mediation[["Pred_Y"]])
    data.table(
      Effect = c("Indirect effect", "Direct effect", "Total effect"),
      Estimate = c(phi_11 - phi_10, phi_10 - phi_00)
    )
  }

  # Check for error
  check_error(data, id_var, base_vars, exposure, time_var, models)

  # Get time length
  time_len <- length(unique(na.omit(data[[time_var]])))

  intervention <- list(always = 1, never = 0, mediation = NULL)

  # Test if mediator model exists
  med_flag <- sapply(models, function(mods) mods$mod_type == "mediator")
  if (!any(med_flag)) {
    stop("Mediator model was not defined.", domain = "causalMed")
  }

  # Run along the models
  if (verbose) {
    cat("\n====== Fitting models =======\n")
  }

  # Run original estimate
  arg_est <- get_args_for(.gformula)
  arg_est$return_fitted <- TRUE
  est_ori <- do.call(.gformula, arg_est)

  # Run bootstrap
  arg_pools <- get_args_for(bootstrap_helper)
  pools <- do.call(bootstrap_helper, arg_pools)

  # Mean value of the outcome at each time point by intervention
  est_out <- data.table::rbindlist(est_ori$gform.data, idcol = "Intervention")
  est_out <- est_out[, list(Est = mean(Pred_Y)), by = c("Intervention")]
  risk_est <- risk_calc(est_ori$gform.data)

  # Get the mean of bootstrap results
  if (R > 1) {
    pools_res <- lapply(pools, function(bt) {
      out <- sapply(bt, function(x) {
        x[, list(Est = mean(Pred_Y))]
      }, simplify = FALSE)
      data.table::rbindlist(out, idcol = "Intervention")
    })
    pools_res <- data.table::rbindlist(pools_res)

    # Calculate Sd and percentile confidence interval
    pools_res <- pools_res[, .(
      Sd = sd(Est), perct_lcl = quantile(Est, 0.025),
      perct_ucl = quantile(Est, 0.975)
    ),
    by = c("Intervention")
    ]

    # Merge all and calculate the normal confidence interval
    est_out <- merge(est_out, pools_res, by = c("Intervention"))
    est_out <- est_out[, `:=`(
      norm_lcl = Est - stats::qnorm(0.975) * Sd,
      norm_ucl = Est + stats::qnorm(0.975) * Sd
    )]

    res_pools <- lapply(pools, risk_calc)
    res_pools <- data.table::rbindlist(res_pools)
    # Calculate Sd and percentile confidence interval
    res_pools <- res_pools[, .(
      Sd = sd(Estimate, na.rm = TRUE),
      perct_lcl = quantile(Estimate, 0.025, na.rm = TRUE),
      perct_ucl = quantile(Estimate, 0.975, na.rm = TRUE)
    ),
    by = c("Effect")
    ]

    # Merge all and calculate the normal confidence interval
    risk_est <- merge(risk_est, res_pools, by = c("Effect"))
    risk_est <- risk_est[, `:=`(
      norm_lcl = Estimate - stats::qnorm(0.975) * Sd,
      norm_ucl = Estimate + stats::qnorm(0.975) * Sd
    )]
  }

  # Extract fitted model information
  resp_vars_list <- sapply(est_ori$fitted.models, function(x) {
    x$rsp_vars
  })
  if (return_fitted) {
    fitted_mods <- lapply(est_ori$fitted.models, function(x) {
      x$fitted
    })
  } else {
    fitted_mods <- lapply(est_ori$fitted.models, function(x) {
      summary(x$fitted)$coefficients
    })
  }
  names(fitted_mods) <- resp_vars_list

  return(list(
    call = tpcall,
    estimate = risk_est,
    risk_size = est_out,
    gform.data = data.table::rbindlist(est_ori$gform.data, idcol = "Intervention"),
    fitted.models = fitted_mods
  ))
}
