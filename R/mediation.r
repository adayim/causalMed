
#' G-formula based mediation analysis
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
                      return_data = FALSE,
                      R = 500) {

  tpcall <- match.call()
  mediation_type <- match.arg(mediation_type)

  # Calculate mediation effect
  risk_calc <- function(data_list, return_data) {
    if(return_data){
      phi_11 <- sum(data_list$always[["Pred_Y"]]) / length(data_list$always[["Pred_Y"]])
      phi_00 <- sum(data_list$never[["Pred_Y"]]) / length(data_list$never[["Pred_Y"]])
      phi_10 <- sum(data_list$mediation[["Pred_Y"]]) / length(data_list$mediation[["Pred_Y"]])
      data.table(Effect = c("Indirect effect", "Direct effect", "Total effect"),
                 Est = c(phi_11 - phi_10, phi_10 - phi_00, phi_11 - phi_00)
                )
    }else{
      data.table(Effect = c("Indirect effect", "Direct effect", "Total effect"),
                 Est = c(data_list$always[["Pred_Y"]] - data_list$mediation[["Pred_Y"]], 
                         data_list$mediation[["Pred_Y"]] - data_list$never[["Pred_Y"]],
                         data_list$always[["Pred_Y"]] - data_list$never[["Pred_Y"]])
                )
    }
    
  }
  risk_calc <- function(data_list, return_data) {
    if(return_data){
      phi_11 <- sum(data_list$always) / length(data_list$always)
      phi_00 <- sum(data_list$never) / length(data_list$never)
      phi_10 <- sum(data_list$mediation) / length(data_list$mediation)
      data.table(Effect = c("Indirect effect", "Direct effect", "Total effect"),
                 Est = c(phi_11 - phi_10, phi_10 - phi_00, phi_11 - phi_00)
      )
    }else{
      data.table(Effect = c("Indirect effect", "Direct effect", "Total effect"),
                 Est = c(data_list$always - data_list$mediation, 
                         data_list$mediation - data_list$never,
                         data_list$always - data_list$never)
      )
    }
    
  }

  if (exists(".Random.seed")) {
    orig.seed <- get(".Random.seed", .GlobalEnv)
    on.exit(.Random.seed <<- orig.seed)
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

  # Run original estimate
  arg_est <- get_args_for(.gformula)
  arg_est$return_fitted <- TRUE
  est_ori <- do.call(.gformula, arg_est)

  # Run bootstrap
  arg_pools <- get_args_for(bootstrap_helper)
  pools <- do.call(bootstrap_helper, arg_pools)

  # Mean value of the outcome at each time point by intervention
  if(return_data){
    est_out <- data.table::rbindlist(est_ori$gform.data, idcol = "Intervention")
    est_out <- est_out[, list(Est = sum(Pred_Y) / length(Pred_Y)), by = c("Intervention")]
  }else{
    est_out <- data.table::as.data.table(utils::stack(est_ori$gform.data))
    colnames(est_out) <- c("Est", "Intervention")
  }
  risk_est <- risk_calc(est_ori$gform.data, return_data = return_data)
  

  # Get the mean of bootstrap results
  if (R > 1) {
    pools_res <- lapply(pools, function(bt) {
      out <- utils::stack(bt)
      colnames(out) <- c("Est", "Intervention")
      return(out)
    })
    pools_res <- data.table::rbindlist(pools_res)

    # Calculate Sd and percentile confidence interval
    pools_res <- pools_res[, .(
      Sd = sd(Est), 
      perct_lcl = quantile(Est, 0.025, na.rm = TRUE),
      perct_ucl = quantile(Est, 0.975, na.rm = TRUE)
    ),
    by = c("Intervention")
    ]

    # Merge all and calculate the normal confidence interval
    est_out <- merge(est_out, pools_res, by = c("Intervention"))
    est_out <- est_out[, `:=`(
      norm_lcl = Est - stats::qnorm(0.975) * Sd,
      norm_ucl = Est + stats::qnorm(0.975) * Sd
    )]

    res_pools <- lapply(pools, risk_calc, return_data = return_data)
    res_pools <- data.table::rbindlist(res_pools)
    # Calculate Sd and percentile confidence interval
    res_pools <- res_pools[, .(
      Sd = sd(Est, na.rm = TRUE),
      perct_lcl = quantile(Est, 0.025, na.rm = TRUE),
      perct_ucl = quantile(Est, 0.975, na.rm = TRUE)
    ),
    by = c("Effect")
    ]

    # Merge all and calculate the normal confidence interval
    risk_est <- merge(risk_est, res_pools, by = c("Effect"))
    risk_est <- risk_est[, `:=`(
      norm_lcl = Est - stats::qnorm(0.975) * Sd,
      norm_ucl = Est + stats::qnorm(0.975) * Sd
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
      list(call = x$fitted$call, coeff = summary(x$fitted)$coefficients)
    })
  }
  names(fitted_mods) <- resp_vars_list

  # Return data
  if(return_data){
    dat_out <- data.table::rbindlist(est_ori$gform.data, idcol = "Intervention")
  }else{
    dat_out <- NULL
  }

  return(list(
    call = tpcall,
    estimate = risk_est,
    effect_size = est_out,
    sim_data = dat_out,
    fitted_models = fitted_mods
  ))
}
