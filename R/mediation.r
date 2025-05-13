
#' G-formula based mediation analysis
#'
#' @description Mediation analysis for time varying mediator, estimation based
#'  on g-formula. Output contains total effect, natural direct effect and natural
#'   indirect effect for mediation or regular g-formula. data.frame will be returned.
#'
#' @note Not that final outcome must be the same for in all rows per subject.
#'   If the dataset is  survival settings, records after the interested outcome event
#'    must be deleted. The function it self do some data manipulation internally.
#'    Please prepare the data as longitudinal (long data) format.
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
                      mc_sample = nrow(data)*100,
                      mediation_type = c("N", "I"),
                      return_fitted = FALSE,
                      return_data = FALSE,
                      R = 500,
                      quiet = FALSE,
                      seed = mc_sample*100) {

  tpcall <- match.call()
  all.args <- mget(names(formals()),sys.frame(sys.nframe()))

  # Initilise warning
  init_warn()

  mediation_type <- match.arg(mediation_type)

  set.seed(seed)

  # Check for error
  check_error(data, id_var, base_vars, exposure, time_var, models)

  # Get time length
  time_len <- length(unique(na.omit(data[[time_var]])))

  intervention <- list(Ph11 = 1, Phi00 = 0, Phi10 = NULL)

  # Test if mediator model exists
  med_flag <- sapply(models, function(mods) mods$mod_type == "mediator")
  if (!any(med_flag)) {
    stop("Mediator model was not defined.", domain = "causalMed")
  }

  # Run original estimate
  arg_est <- get_args_for(.gformula)
  arg_est$return_fitted <- TRUE
  arg_est$progress_bar <- substitute(quiet, env = parent.frame())
  est_ori <- do.call(.gformula, arg_est)

  # Mean value of the outcome at each time point by intervention
  if(return_data){
    est_out <- data.table::rbindlist(est_ori$gform.data, idcol = "Intervention")
    est_out <- est_out[, list(Est = sum(Pred_Y) / length(Pred_Y)), by = c("Intervention")]
  }else{
    est_out <- data.table::as.data.table(utils::stack(est_ori$gform.data))
    colnames(est_out) <- c("Est", "Intervention")
  }

  if(return_data)
    estimate_extract <- sapply(est_ori$gform.data, "[[", "Pred_Y", simplify = FALSE)
  else
    estimate_extract <- est_ori$gform.data

  risk_est <- risk_estimate.mediation(estimate_extract, return_data = return_data)


  # Get the mean of bootstrap results
  if (R > 1) {
    # Run bootstrap
    arg_pools <- get_args_for(bootstrap_helper)
    pools <- do.call(bootstrap_helper, arg_pools)

    pools_res <- lapply(pools, function(bt) {
      out <- utils::stack(bt$gform.data)
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
    est_out <- merge(est_out, pools_res, by = c("Intervention"), sort = FALSE)
    est_out <- est_out[, `:=`(
      norm_lcl = Est - stats::qnorm(0.975) * Sd,
      norm_ucl = Est + stats::qnorm(0.975) * Sd
    )]

    res_pools <- lapply(pools, function(x) x$gform.data)
    res_pools <- lapply(res_pools, risk_estimate.mediation, return_data = return_data)
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
    risk_est <- merge(risk_est, res_pools, by = c("Effect"), sort = FALSE)
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

  if(length(causalmed_env$warning)>0){
    message(paste(causalmed_env$warning, collapse = "\n=============\n"),
            domain = "causalMed")
  }

  y <- list(call = tpcall,
            all.args = all.args,
            estimate = risk_est,
            effect_size = est_out,
            sim_data = dat_out,
            fitted_models = fitted_mods
          )
  class(y) <- c("gformula", class(y))
  return(y)
}
